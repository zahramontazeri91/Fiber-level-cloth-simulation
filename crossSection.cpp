#include "CrossSection.h"
#include "curveFitting.h"
#include <Eigen/SVD>

#if 0
/* constructor for procedural yarn (for debug use) */
CrossSection::CrossSection(const Fiber::Yarn &yarn) {
	m_yarn = yarn;
	m_planesList.resize(yarn.getStepNum());
	for (int step_id = 0; step_id < m_planesList.size(); step_id++) {
		const float z = yarn.getStepSize() * (step_id - yarn.getStepNum() / 2.f); // devided by 2 Bcuz yarn lies between neg and pos z
		m_planesList[step_id].point = vec3f(0.f, 0.f, z);
		m_planesList[step_id].n = vec3f(0.f, 0.f, 1.f); //always parallel to xy plane
		m_planesList[step_id].e1 = vec3f(1.f, 0.f, 0.f);
		m_planesList[step_id].e2 = vec3f(0.f, 1.f, 0.f);
	}
	// no need to initialize m_curve since we already have the cross-section planes
}
#endif

void CrossSection::init(const char* yarnfile, const int ply_num, const char* curvefile,
	const int seg_subdiv, const int num_planes, std::vector<yarnIntersect2D> &allPlaneIntersect) {
	m_yarn.build(yarnfile, ply_num);
	//TODO: pass this to the func
	m_curve.init_norm(curvefile,"write2normals.txt", seg_subdiv);
	std::vector<yarnIntersect> itsLists;
	buildPlanes(num_planes, itsLists);
	std::cout << "Finding the intersections with the cross-sectional plane... \n";
	PlanesIntersections2D(itsLists, allPlaneIntersect);
}
void CrossSection::init(const char* yarnfile, const int ply_num, const char* curvefile, const char* normfile,
	const int seg_subdiv, const int num_planes, std::vector<yarnIntersect2D> &allPlaneIntersect) {
	m_yarn.build(yarnfile, ply_num);
	m_curve.init(curvefile, normfile, seg_subdiv);
	std::vector<yarnIntersect> itsLists;
	buildPlanes(num_planes, itsLists);
	std::cout << "Finding the intersections with the cross-sectional plane... \n";
	PlanesIntersections2D(itsLists, allPlaneIntersect);
}

void CrossSection::buildPlanes(const int num_planes, std::vector<yarnIntersect> &itsLists) {
	const double curveLen = m_curve.totalLength();
	//const double crossSectionLen = curveLen / static_cast<double>(num_planes - 1); //place plane at the very ends as well
	const double crossSectionLen = curveLen / static_cast<double>(num_planes + 1); //don't place planes at the end-points
	const double crossSection_t = m_curve.arcLengthInvApprox(crossSectionLen);
	m_planesList.resize(num_planes);
	for (int i = 0; i < num_planes; ++i) {
		//float current_t = i * crossSection_t;//place plane at the very ends as well
		float current_t = (i + 1) * crossSection_t;
		Eigen::Vector3d curve_p = m_curve.eval(current_t);
		Eigen::Vector3d curve_t = m_curve.evalTangent(current_t);
		/* need binormal and normal of the plane to project 3D points on the plane */
		Eigen::Vector3d curve_n = m_curve.evalNormal(current_t); //short axis
		Eigen::Vector3d curve_b = curve_t.cross(curve_n); //long axis
		m_planesList[i].point = vec3f(curve_p[0], curve_p[1], curve_p[2]);
		m_planesList[i].n = vec3f(curve_t[0], curve_t[1], curve_t[2]);
		m_planesList[i].e1    = vec3f(curve_n[0], curve_n[1], curve_n[2]);
		m_planesList[i].e2 = cross(m_planesList[i].n, m_planesList[i].e1);
		//assert(dot(m_planesList[i].n, m_planesList[i].e1) < EPS  && "n and e1 are not perpendicular\n");
	}
	std::vector<std::vector<vec3f>> plyCenters;
	allPlyCenters(plyCenters, itsLists);
	// Define inplane 2D coord using the direction from yarn-center to intersection of first ply-center 
	for (int i = 0; i < num_planes; ++i) {
		// plyCenters[i][0] - m_planesList[i].point is too small so its dot product with n doesn't show they are perpendicular (e1.n !=0 )
		//m_planesList[i].e1 = nv::normalize(plyCenters[i][0] - m_planesList[i].point);
		//m_planesList[i].e2 = cross(m_planesList[i].n, m_planesList[i].e1);
	}

	//debug: testing the e# coord
	//std::vector<std::vector<vec3f>> plyCenters;
	//allPlyCenters(plyCenters);
	//// Define inplane 2D coord using the direction from yarn-center to intersection of first ply-center 
	//FILE *fout;
	//if (fopen_s(&fout, "../data/test_planeCoord.txt", "wt") == 0) {
	//	for (int i = 0; i < num_planes; ++i) {
	//		if (i % 10 == 0) {
	//			fprintf_s(fout, "%.4f %.4f %.4f ", m_planesList[i].point.x, m_planesList[i].point.y, m_planesList[i].point.z);
	//			fprintf_s(fout, "%.4f %.4f %.4f ", m_planesList[i].n.x, m_planesList[i].n.y, m_planesList[i].n.z);
	//			fprintf_s(fout, "%.4f %.4f %.4f ", m_planesList[i].e1.x, m_planesList[i].e1.y, m_planesList[i].e1.z);
	//			fprintf_s(fout, "\n");
	//		}
	//	}
	//	fclose(fout);
	//} 
}

bool CrossSection::linePlaneIntersection(const vec3f &start, const vec3f &end, const Plane &plane, vec3f &its) {
	const vec3f dir = nv::normalize(end - start);
	const float t_end = length(end - start);
	const float denom = dot(plane.n, dir);
	if (std::abs(denom) > EPS) {
		vec3f cntrStart = plane.point - start;
		float t = dot(cntrStart, plane.n) / denom;
		if (t >= 0.f && t <= t_end) {
			its = start + dir*t;

			assert(length(its - plane.point) > 1e-6 && "intersection exactly at the plane.point"); //TODO: handle this corner case later 
			assert(std::abs(dot(its - plane.point, plane.n)) < 1e-6);
			return true;
		}
	}
	return false;
}

bool CrossSection::yarnPlaneIntersection(const Plane &plane, yarnIntersect &itsList) {
	bool isIntrsct = false;
	const int ply_num = m_yarn.plys.size();
	itsList.resize(ply_num);

	assert(nv::length(plane.n) - 1.f < EPS   && "normal vector is not normalized!");

	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)

	for (int p = 0; p < ply_num; ++p) {
		const int fiber_num = m_yarn.plys[p].fibers.size();
		for (int f = 0; f < fiber_num; ++f) {
			int hitFiberNum = 0;
			const int vrtx_num = m_yarn.plys[p].fibers[f].vertices.size();
			float min_dist = std::numeric_limits<float>::max();
			vec3f closest_its(0.f);
			for (int v = 0; v < vrtx_num - 1; ++v) { //Each segment is v[i] to v[i+1]
				vec3f start = m_yarn.plys[p].fibers[f].vertices[v];
				vec3f end = (m_yarn.plys[p].fibers[f].vertices[v + 1]);
				vec3f its(0.f);
				if (linePlaneIntersection(start, end, plane, its)) {

					if (hitFiberNum) {
						hitFiberNum++;
						float dist = nv::distance(its, plane.point); //distance between the fiber and the center
						if (dist < min_dist)
							closest_its = its;
						//std::cout << "Ply " << p << ", fiber " << f << " intersects " << hitFiberNum + 1 << " times with the plane! \n";
					}
					isIntrsct = true;
					closest_its = its;
					hitFiberNum++;
				}
			}
			if (hitFiberNum) 
				itsList[p].push_back(closest_its);
		}
	}

	if (isIntrsct)
		return true;
	return false;
}

bool CrossSection::allPlanesIntersections(std::vector<yarnIntersect> &itsLists) {
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores) 

	bool isIntrsct = false;
	for (int i = 0; i < m_planesList.size(); ++i) {
		yarnIntersect itsList;
		if (yarnPlaneIntersection(m_planesList[i], itsList)) {
			isIntrsct = true;
			itsLists.push_back(itsList);
		}
	}

	if (isIntrsct)
		return true;
	return false;
}

void CrossSection::shapeMatch(const Eigen::MatrixXf &pnt_trans, const Eigen::MatrixXf &pnt_ref, Ellipse &ellipse, float &theta_R) {
	
	assert(pnt_trans.cols() == pnt_ref.cols());
	const int n = pnt_trans.cols();
	Eigen::Matrix2f Apq = Eigen::Matrix2f::Zero();
	Eigen::Matrix2f Aqq_1 = Eigen::Matrix2f::Zero();
	Eigen::Matrix2f Aqq = Eigen::Matrix2f::Zero();

	Eigen::MatrixXf Qtrans = pnt_ref.transpose();
	Eigen::MatrixXf cm_ref = Eigen::MatrixXf::Zero(2,1);
	Eigen::MatrixXf cm_trans = Eigen::MatrixXf::Zero(2,1);
	if (1) { //allow translation
		for (int i = 0; i < n; ++i) {
			cm_ref += pnt_ref.col(i);
			cm_trans += pnt_trans.col(i);
		}
		cm_ref /= n;
		cm_trans /= n;
	}

	for (int i = 0; i < n; ++i) {
		Apq	  += (pnt_trans.col(i) - cm_trans)  *  (pnt_ref.col(i) - cm_ref).transpose() / n;
		Aqq_1 += (pnt_ref.col(i) - cm_ref) *  (pnt_ref.col(i) - cm_ref).transpose() / n;
	}

	assert(Aqq_1.determinant() != 0);
	Aqq = Aqq_1.inverse();

	//SVD decomposition
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(Apq*Aqq, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::Matrix2f U = svd.matrixU();
	Eigen::Matrix2f sigma = svd.singularValues().asDiagonal();
	Eigen::Matrix2f V = svd.matrixV();

	//RS decompose
	Eigen::Matrix2f S = V*sigma*V.transpose();
	Eigen::Matrix2f R = U*V.transpose();
	ellipse.longR = std::max(sigma(0,0), sigma(1,1));
	ellipse.shortR = std::min(sigma(0, 0), sigma(1, 1));
	//Make V reflection free (if det(V) is negative)
	Eigen::Matrix2f m;
	m << 1, 0, 0, -1;
	if (V.determinant() < 0)
		V = V*m;
	ellipse.angle = atan2(V(1, 0), V(0, 0));
	ellipse.angle = ellipse.angle < 0 ? ellipse.angle + 2.f*pi : ellipse.angle;

	theta_R = atan2(R(1, 0), R(0, 0));
	theta_R = theta_R < 0 ? theta_R + 2.f*pi : theta_R;
}

void CrossSection::yarnShapeMatch(const yarnIntersect2D &pnts_trans, const yarnIntersect2D &pnts_ref, Ellipse &ellipse, float &theta_R) {

	//Find the total number of points for all plys
	int sz_ref = 0, sz_trans = 0;
	for (int p = 0; p < pnts_ref.size(); ++p)
		sz_ref += pnts_ref[p].size();
	for (int p = 0; p < pnts_trans.size(); ++p)
		sz_trans += pnts_trans[p].size();

	assert(sz_ref == sz_trans);
	const int n = sz_ref;
	Eigen::MatrixXf all_ref(2, n);
	Eigen::MatrixXf all_trans(2, n);
	//if (sz_ref != sz_trans) {
	//	all_ref = Eigen::MatrixXf::Random(2, n);
	//	all_trans = Eigen::MatrixXf::Random(2, n);
	//	shapeMatch(all_trans, all_ref, ellipse, theta_R);
	//	return;
	//}

	int c = 0;
	for (int p = 0; p < pnts_ref.size(); ++p) {
		for (int i = 0; i < pnts_ref[p].size(); ++i) {
			all_ref(0, c) = pnts_ref[p][i].x;
			all_ref(1, c) = pnts_ref[p][i].y;
			++c;
		}
	}


	c = 0;
	for (int p = 0; p < pnts_trans.size(); ++p) {
		for (int i = 0; i < pnts_trans[p].size(); ++i) {
			all_trans(0, c) = pnts_trans[p][i].x;
			all_trans(1, c) = pnts_trans[p][i].y;
			++c;
		}
	}

	shapeMatch(all_trans, all_ref, ellipse, theta_R);	
}
void CrossSection::yarnShapeMatches(const std::vector<yarnIntersect2D> &pnts_trans, const std::vector<yarnIntersect2D> &pnts_ref, std::vector<Ellipse> &ellipses, std::vector<float> &all_theta_R) {
	assert(pnts_trans.size() == pnts_trans.size());
	const int n = pnts_trans.size();
	ellipses.resize(n);
	all_theta_R.resize(n);

	for (int i = 0; i < n; ++i) {
		Ellipse ellipse;
		float theta_R;
		if (pnts_trans[i][0].size() != pnts_ref[i][0].size() || pnts_trans[i][1].size() != pnts_ref[i][1].size()) {
			//std::cout << i << " is not valid cross-section"  <<  "\n";
			//std::cout << pnts_trans[i][0].size() << " " << pnts_ref[i][0].size() << std::endl;
			//std::cout << pnts_ref[i][0].size() << " " << pnts_ref[i][0].size() << std::endl;
			ellipses[i].shortR = 0.0;
			ellipses[i].longR = 0.0;
			all_theta_R[i] = all_theta_R[i-1];
			continue;
		}

		yarnShapeMatch(pnts_trans[i], pnts_ref[i], ellipse, theta_R);
		ellipses[i] = ellipse;
		all_theta_R[i] = theta_R;

	}
}

void CrossSection::allPlyCenters(std::vector<std::vector<vec3f>> &plyCenters, std::vector<yarnIntersect> &itsLists) {

	allPlanesIntersections(itsLists);

	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)

	const int plane_num = itsLists.size();
	plyCenters.resize(plane_num); //number of planes
	for (int i = 0; i < plane_num; ++i) { //plane num
		const int ply_num = itsLists[i].size();
		plyCenters[i].resize(ply_num);
		for (int p = 0; p < ply_num; ++p) { //ply num
			const int fiber_num = itsLists[i][p].size();
			for (int v = 0; v < fiber_num; ++v) { // fiber-centers num
				plyCenters[i][p] += itsLists[i][p][v];
			}
			plyCenters[i][p] /= fiber_num;
		}
	}

	//for testing:
	/*FILE *fout;
	if (fopen_s(&fout, "../data/test_plyCenter.txt", "wt") == 0) {
	fprintf_s(fout, "%d\n", plane_num);
	for (int i = 0; i < plane_num; ++i) {
	fprintf_s(fout, "%.4f %.4f %.4f \n", plyCenters[i][0].x, plyCenters[i][0].y, plyCenters[i][0].z);
	}
	fclose(fout);
	}*/
}

void CrossSection::write_PlanesIntersections3D(const char* filename, std::vector<yarnIntersect> &itsLists) {

	assert(m_planesList.size() == itsLists.size());
	//std::cout << itsLists.size() << " out of " << m_planesList.size()  << " many planes had intersections! \n";
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)

	const int ply_num = m_yarn.plys.size();
	FILE *fout;
	if (fopen_s(&fout, filename, "wt") == 0) {

		fprintf_s(fout, "plane_num: %d \n", itsLists.size());
		fprintf_s(fout, "ply_num: %d \n", ply_num);
		fprintf_s(fout, "\n");

		for (int cs = 0; cs < itsLists.size(); ++cs) { 
			/* write all the intersections for each ply */
			for (int p = 0; p < ply_num; ++p) { //for each ply 
				fprintf_s(fout, "ply_fiber_num: %d \n", itsLists[cs][p].size());
				fprintf_s(fout, "plyCenter: %.4lf %.4lf %.4lf \n", itsLists[cs][p][0].x, itsLists[cs][p][0].y, itsLists[cs][p][0].z);
				for (int j = 0; j < itsLists[cs][p].size(); ++j) { //for each intersected fiber
					vec3f fiberIts = itsLists[cs][p][j];
					fprintf_s(fout, "%.4lf %.4lf %.4lf \n", fiberIts.x, fiberIts.y, fiberIts.z);
				}
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
	//std::cout << "Intersections are written to the file successfully! \n\n";
}

void CrossSection::project2Plane(const vec3f& P3d, const Plane& plane, vec2f& P2d) {
	vec3f n = plane.n;
	vec3f e1 = plane.e1;
	vec3f e2 = plane.e2;

	assert(length(e1) - 1.f < 1e-6);
	assert(length(e2) - 1.f < 1e-6);
	assert(dot(e1, e2) < EPS);
	P2d.x = dot(e1, (P3d - plane.point));
	P2d.y = dot(e2, (P3d - plane.point));
}

void CrossSection::PlanesIntersections2D(std::vector<yarnIntersect> &itsLists, std::vector<yarnIntersect2D> &allPlaneIntersect) {
	if (m_planesList.size() != itsLists.size())
		std::cout << itsLists.size() << " out of " << m_planesList.size() << " many planes had intersections! \n";
	assert(m_planesList.size() == itsLists.size() && "not all the planes have intersections!");

	const int ply_num = m_yarn.plys.size();
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)

	for (int cs = 0; cs < itsLists.size(); ++cs) { //for each plane
		Plane plane;
		get_plane(cs, plane);
		/* First write the yarn-center */
		vec2f center2D;
		project2Plane(plane.point, plane, center2D);
		assert(center2D == vec2f(0.f));

		/* Then all the intersections for each ply */
		yarnIntersect2D planeProjected;
		planeProjected.resize(ply_num);
		for (int p = 0; p < ply_num; ++p) { //for each ply 
			for (int j = 0; j < itsLists[cs][p].size(); ++j) { //for each intersected fiber
				vec3f fiberIts = itsLists[cs][p][j];
				vec2f projected;
				project2Plane(fiberIts, plane, projected);
				planeProjected[p].push_back(projected);
			}
		}
		allPlaneIntersect.push_back(planeProjected);
	}
	//std::cout << "Intersections are written to the file successfully! \n\n";
}


void findEigenValue(const cv::Mat &data_pts, const std::vector<vec2f> &eigen_vecs, std::vector<float> &eigen_val) {
	float max0 = std::numeric_limits<float>::min();
	float max1 = std::numeric_limits<float>::min();
	for (int i = 0; i < data_pts.rows; ++i) {
		float prj0 = eigen_vecs[0].x * data_pts.at<float>(i, 0) + eigen_vecs[0].y * data_pts.at<float>(i, 1);
		if (std::abs(prj0) > max0)
			max0 = std::abs(prj0);
		float prj1 = eigen_vecs[1].x * data_pts.at<float>(i, 0) + eigen_vecs[1].y * data_pts.at<float>(i, 1);
		if (std::abs(prj1) > max1)
			max1 = std::abs(prj1);
	}

	eigen_val[0] = max0;
	eigen_val[1] = max1;
}

void CrossSection::preComputeEllipses(const std::vector<Ellipse> &ellipses, Eigen::MatrixXf &R1, Eigen::MatrixXf &R2, Eigen::MatrixXf &theta) {
	const int plane_num = ellipses.size();
	R1.resize(plane_num, 4);
	R2.resize(plane_num, 4);
	theta.resize(plane_num, 4);

	for (int cs = 0; cs < plane_num; ++cs) {

		//assign to R1
		for (int i = 0; i < 4; ++i) {
			R1(cs, i) = i % 2 ? ellipses[cs].shortR : ellipses[cs].longR;
		}
		//assign to R2
		for (int i = 0; i < 4; ++i) {
			R2(cs, i) = i % 2  ? ellipses[cs].longR : ellipses[cs].shortR;
		}
		//assign to theta
		float angle = ellipses[cs].angle;
		theta(cs, 0) = angle;
		for (int i = 1; i < 4; ++i) {
			angle += pi / 2.f;
			angle = angle > 2.f*pi ? angle - 2.f*pi : angle;
			theta(cs, i) = angle;
		}	
	}
}


void CrossSection::costFunction(const Eigen::MatrixXf &R1, const Eigen::MatrixXf &R2, const Eigen::MatrixXf &theta, const std::vector<bool> &isValid, std::vector<Eigen::Matrix4f> &cost) {
	assert(isValid[0]); //assumption
	Eigen::Matrix4f config = Eigen::Matrix4f::Zero();
	// for the ith cross-section, last valid one is known using isValid()
	// so config is defined as a 4x4 matrix for all configs for two consecutive planes
	cost.resize(isValid.size());
	cost[0] = config;
	for (int i = 1; i < cost.size(); ++i) { //started from 1
		cost[i] = Eigen::Matrix4f::Zero(); //keep it zero for in-valid ones
		if (isValid[i]) {
			int ip = i - 1;
			while (!isValid[ip])
				ip = ip - 1; //ip is the last valid cs
			for (int k_i = 0; k_i < 4; ++k_i) {
				for (int k_ip = 0; k_ip < 4; ++k_ip) {
					float w1 = 1.f / pow(0.05,2); //yarn radius
					float w2 = 1.f / pow(0.05,2);
					float w3 = 1.f / pow(2.f*pi, 2);
					float dR1 = (R1(i, k_i) - R1(ip, k_ip)) / (i - ip);
					float dR2 = (R2(i, k_i) - R2(ip, k_ip)) / (i - ip);
					float dtheta = std::abs(theta(i, k_i) - theta(ip, k_ip) );
					dtheta = std::min(dtheta, 2.f*pi - dtheta)/(i - ip);
					float g = (i - ip) * ( w1*pow(dR1, 2) + w2*pow(dR1,2) + w3*pow(dtheta, 2) );
					config(k_i, k_ip) = g;
				}
			}	
			cost[i] = config;
		}	
	}
}

void CrossSection::dynamicProgramming(const std::vector<bool> &isValid, const std::vector<Eigen::Matrix4f> &cost, Eigen::MatrixXf &totalCost, Eigen::MatrixXf &preConfig) {
	const int plane_num = isValid.size();
	totalCost.resize(plane_num, 4);
	preConfig.resize(plane_num, 4);
	//initialization for first cs
	for (int i = 0; i < 4; ++i) {
		totalCost(0, i) = 0.f;
		preConfig(i, i) = -1.f;
	}

	for (int i = 1; i < plane_num; ++i) { //started from 1
		if (isValid[i]) {
			int ip = i - 1;
			while (!isValid[ip])
				ip = ip - 1; //ip is the last valid cs
			for (int k_i = 0; k_i < 4; ++k_i) {
				totalCost(i,k_i) = std::numeric_limits<float>::max();
				for (int k_ip = 0; k_ip < 4; ++k_ip) {
					float val = totalCost(ip, k_ip) + cost[i](k_i, k_ip);
					if (val < totalCost(i, k_i)) {
						totalCost(i, k_i) = val;
						preConfig(i, k_i) = k_ip;
					}
				}
			}
		}
		//totalCost and preconfig for in-valid ones are not assigned
	}
}

void interpolateEllipses(const std::vector<Ellipse> &ellipses, const std::vector<bool> &isValid, std::vector<Ellipse> &ellipses_smooth) {
	ellipses_smooth.resize(ellipses.size() );
	assert(ellipses.size() == isValid.size() );
	assert(isValid[0]); //assumption: handle this later
	int i = 0;
	while (i < isValid.size()) {
		if (!isValid[i]) { // i-1 is valid 
			int j = i+1;
			while (!isValid[j]) // j is valid
				j++;
			//interpolate values between i-1 and j 
			float delta_indx = (j) - (i - 1) - 1;
			float delta_longR = ellipses[j ].longR - ellipses[i - 1].longR;
			float delta_shortR = ellipses[j ].shortR - ellipses[i - 1].shortR;
			float delta_angle = ellipses[j ].angle - ellipses[i - 1].angle;
			int cnt = 1;
			for (int cs = i; cs < j; ++cs) {
				ellipses_smooth[cs].longR = ellipses[i - 1].longR + cnt * delta_longR / delta_indx;
				ellipses_smooth[cs].shortR = ellipses[i - 1].shortR + cnt * delta_shortR / delta_indx;
				ellipses_smooth[cs].angle = ellipses[i - 1].angle + cnt * delta_angle / delta_indx;
				cnt++;
			}
			//interpolate from other direction
			if (std::abs(delta_angle) > 2.f*pi - std::abs(delta_angle)) { 
				delta_angle = 2.f*pi - std::abs(ellipses[j].angle - ellipses[i - 1].angle);
				delta_angle = ellipses[j].angle > ellipses[i - 1].angle ? -delta_angle : delta_angle;
				int cnt = 1;
				for (int cs = i; cs < j; ++cs) {
					ellipses_smooth[cs].angle = ellipses[i - 1].angle + cnt * delta_angle / delta_indx;
					if (ellipses_smooth[cs].angle > 2 * pi)
						ellipses_smooth[cs].angle -= 2 * pi;
					if (ellipses_smooth[cs].angle < 0)
						ellipses_smooth[cs].angle += 2 * pi;
					cnt++;
				}
			}
			i = j;
		}
		else {
			ellipses_smooth[i].longR = ellipses[i].longR;
			ellipses_smooth[i].shortR = ellipses[i].shortR;
			ellipses_smooth[i].angle = ellipses[i].angle;
			++i;
		}
	}
}

void CrossSection::fitEllipses(const std::vector<yarnIntersect2D> &allpts, std::vector<Ellipse> &ellipses, std::vector<bool> &isValid) {

	const int plane_num = allpts.size();
	isValid.resize(plane_num, true);
	ellipses.resize(plane_num);

	for (int cs = 0; cs < plane_num; ++cs) {
		Ellipse ellipse;
		//Find the total number of points for all plys
		int sz = 0;
		yarnIntersect2D pts = allpts[cs];
		for (int p = 0; p < pts.size(); ++p)
			sz += pts[p].size();
		cv::Mat data_pts(sz, 2, CV_32FC1, cv::Scalar::all(0));

		int c = data_pts.rows;
		for (int p = 0; p < pts.size(); ++p) {
			for (int i = 0; i < pts[p].size(); ++i) {
				data_pts.at<float>(c, 0) = pts[p][i].x;
				data_pts.at<float>(c, 1) = pts[p][i].y;
				--c;
			}
		}

		//Perform PCA analysis
		cv::PCA pca_analysis(data_pts, cv::Mat(), cv::PCA::DATA_AS_ROW, 2);

		//Store the center of the object
		ellipse.center = vec2f(pca_analysis.mean.at<float>(0, 0), pca_analysis.mean.at<float>(0, 1));

		//Store the eigenvalues and eigenvectors
		std::vector<vec2f> eigen_vecs(2);
		std::vector<float> eigen_val(2);
		for (int i = 0; i < 2; ++i)
		{
			eigen_vecs[i] = vec2f(pca_analysis.eigenvectors.at<float>(i, 0),
				pca_analysis.eigenvectors.at<float>(i, 1));
			//eigen_val[i] = pca_analysis.eigenvalues.at<float>(0, i); // Wrong values
		}
		assert(!dot(eigen_vecs[0], eigen_vecs[1]) && "Eigen vectors aren't orthogonal!");

		//find eigen values by projecting the points on eigenVectors
		findEigenValue(data_pts, eigen_vecs, eigen_val);

		//assign to R
		ellipse.longR = eigen_val[0];
		ellipse.shortR = eigen_val[1];


		//assign to theta
		float angle = atan2(eigen_vecs[0].y, eigen_vecs[0].x);
		ellipse.angle = angle < 0 ? angle + 2.f*pi : angle;

		ellipses[cs] = ellipse;

		//assign to isValid
		if (std::abs(eigen_val[0] - eigen_val[1]) < std::max(eigen_val[0], eigen_val[1]) / 2.f) {
			isValid[cs] = false;
		}
		
	}
	std::cout << "here in fitellipse " << ellipses[10].longR << std::endl;
}

#if 0
void CrossSection::fitEllipse(const yarnIntersect2D &pts, Ellipse &ellipse, vec2f &axis1_old, vec2f &axis1_new, const int plane_indx)
{
	//Find the total number of points for all plys
	int sz = 0;
	for (int p = 0; p < pts.size(); ++p)
		sz += pts[p].size();
	cv::Mat data_pts(sz, 2, CV_32FC1, cv::Scalar::all(0));

	int c = data_pts.rows;
	for (int p = 0; p < pts.size(); ++p) {
		for (int i = 0; i < pts[p].size(); ++i) {
			data_pts.at<float>(c, 0) = pts[p][i].x;
			data_pts.at<float>(c, 1) = pts[p][i].y;
			--c;
		}
	}

	//Perform PCA analysis
	cv::PCA pca_analysis(data_pts, cv::Mat(), cv::PCA::DATA_AS_ROW, 2);

	//Store the center of the object
	ellipse.center = vec2f(pca_analysis.mean.at<float>(0, 0), pca_analysis.mean.at<float>(0, 1));


	//Store the eigenvalues and eigenvectors
	std::vector<vec2f> eigen_vecs(2);
	std::vector<float> eigen_val(2);
	for (int i = 0; i < 2; ++i)
	{
		eigen_vecs[i] = vec2f(pca_analysis.eigenvectors.at<float>(i, 0),
			pca_analysis.eigenvectors.at<float>(i, 1));
		//eigen_val[i] = pca_analysis.eigenvalues.at<float>(0, i); // Wrong values
	}

	assert(!dot(eigen_vecs[0], eigen_vecs[1]) && "Eigen vectors aren't orthogonal!");

	//check which eigenVec is axis1 (using minimum rotation algo)
	vec2f axis2_new;

	if (!plane_indx) //first plane
	{
		axis1_old = eigen_vecs[0];
		axis1_new = eigen_vecs[0];
		axis2_new = eigen_vecs[1];
	}
	else {
		//compare the angle between two new axis with the axis1_old

		float angle0 = std::min( pi - acos (dot(eigen_vecs[0], axis1_old)) , acos(dot(eigen_vecs[0], axis1_old)) );
		float angle1 = std::min( pi - acos (dot(eigen_vecs[1], axis1_old)) , acos(dot(eigen_vecs[1], axis1_old)) );

		if (dot(eigen_vecs[0], axis1_old) >= 1 || dot(eigen_vecs[1], axis1_old) >= 1)
		{//handle nan (when axis-new and old are the same)
			axis1_new = axis1_old;
			axis2_new = vec2f(-1 * axis1_new.y, axis1_new.x);
		}
		else if (angle0 < angle1) {
			if (pi - acos(dot(eigen_vecs[0], axis1_old)) < acos(dot(eigen_vecs[0], axis1_old))) //flipped 
			{
				eigen_vecs[0] *= -1;
			}
			axis1_new = eigen_vecs[0];
			axis2_new = vec2f(-1 * axis1_new.y, axis1_new.x);

		}
		else {
			if (pi - acos(dot(eigen_vecs[1], axis1_old)) < acos(dot(eigen_vecs[1], axis1_old))) //flipped
			{
				eigen_vecs[1] *= -1;
			}
			axis1_new = eigen_vecs[1];
			axis2_new = vec2f(-1 * axis1_new.y, axis1_new.x);
		}
	}

	//find eigen values by projecting the points on eigenVectors
	float max0 = std::numeric_limits<float>::min();
	float max1 = std::numeric_limits<float>::min();
	for (int i = 0; i < data_pts.rows; ++i) {
		float prj0 = axis1_new.x * data_pts.at<float>(i, 0) + axis1_new.y * data_pts.at<float>(i, 1);
		if (std::abs(prj0) > max0)
			max0 = std::abs(prj0);
		float prj1 = axis2_new.x * data_pts.at<float>(i, 0) + axis2_new.y * data_pts.at<float>(i, 1);
		if (std::abs(prj1) > max1)
			max1 = std::abs(prj1);
	}
	// TODO: find the eigen values using eigen vectors
	eigen_val[0] = max0;
	eigen_val[1] = max1;
	ellipse.longR = eigen_val[0];
	ellipse.shortR = eigen_val[1];
	ellipse.angle = atan2(axis1_new.y, axis1_new.x); // orientation in radians

	//vec2f p1 = ellipse.center + vec2f(axis1_new.x * eigen_val[0], axis1_new.y * eigen_val[0]);
	//vec2f p2 = ellipse.center - vec2f(axis2_new.x * eigen_val[1], axis2_new.y * eigen_val[1]);
	//vec2f long_axis = p1 - ellipse.center;
	//vec2f short_axis = p2 - ellipse.center;
}
#endif

#if 0
/* Given a ellipse, find the minimum area ellipse that covers 95% of fiber centers (search around the given ellipse) */
void CrossSection::minAreaEllipse(const yarnIntersect2D &pts, const Ellipse &ellipse, Ellipse &minEllipse) {
	//Find the total number of points for all plys
	int sz = 0;
	for (int p = 0; p < pts.size(); ++p)
		sz += pts[p].size();
	std::vector<vec2f> data_pts(sz);

	int c = 0;
	for (int p = 0; p < pts.size(); ++p) {
		for (int i = 0; i < pts[p].size(); ++i) {
			data_pts[c].x = pts[p][i].x;
			data_pts[c].y = pts[p][i].y;
			++c;
		}
	}
	//check the percentage of covered dataPoints

	float Rx = ellipse.longR;
	float Ry = ellipse.shortR;
	//float threshold = 0.001;  //// TODO
	//if (!ellipse.shortR) 
	//	Rx = threshold;
	//if (!ellipse.longR) 
	//	Ry = threshold;

	vec2f cnt = ellipse.center;
	const float stepX = Rx * 0.01;
	const float stepY = Ry * 0.01;

	float minArea = std::numeric_limits<float>::max();
	int kk = 0;
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	// Search for all ellipses between R/2 to 3R/3 to speed up
	float rx = 0.5 * Rx;
	while (rx <= 1.5 * Rx) {
		float ry = 0.5 * Ry;
		while (ry <= 1.5 * Ry) {
			kk++;
			int insidePnt = 0;
			for (int i = 0; i < sz; ++i) {
				float x = data_pts[i].x;
				float y = data_pts[i].y;
				float t = std::pow(std::cos(ellipse.angle)*(x - cnt.x) + std::sin(ellipse.angle)*(y - cnt.y), 2) / std::pow(rx, 2);
				t += std::pow(std::sin(ellipse.angle)*(x - cnt.x) - std::cos(ellipse.angle)*(y - cnt.y), 2) / std::pow(ry, 2);
				if (t <= 1)
					insidePnt++;
			}
			float percent = static_cast<float> (insidePnt) / static_cast<float> (sz);
			float area = rx * ry * pi;

			if (percent > 0.95 && area < minArea) {
				minArea = area;
				minEllipse.longR = rx;
				minEllipse.shortR = ry;
			}
			ry += stepY;
		} // end for ry
		rx += stepX;
	} //end for rx

	minEllipse.center = ellipse.center;
	minEllipse.angle = ellipse.angle;

	if (!minEllipse.shortR)
		minEllipse.shortR = Ry;
	if (!minEllipse.longR)
		minEllipse.longR = Rx;
	assert(minEllipse.shortR && "Ellipse length is zero");
	assert(minEllipse.longR && "Ellipse length is zero");
}
#endif

void preComputeEllipses(const std::vector<yarnIntersect2D> &allpts, Eigen::MatrixXf &R1, Eigen::MatrixXf &R2, Eigen::MatrixXf &theta, std::vector<bool> &isValid) {
	const int plane_num = allpts.size();
	R1.resize(plane_num, 4);
	R2.resize(plane_num, 4);
	theta.resize(plane_num, 4);
	isValid.resize(plane_num, true);

	for (int cs = 0; cs < plane_num; ++cs) {
		Ellipse ellipse;
		//Find the total number of points for all plys
		int sz = 0;
		yarnIntersect2D pts = allpts[cs];
		for (int p = 0; p < pts.size(); ++p)
			sz += pts[p].size();
		cv::Mat data_pts(sz, 2, CV_32FC1, cv::Scalar::all(0));

		int c = data_pts.rows;
		for (int p = 0; p < pts.size(); ++p) {
			for (int i = 0; i < pts[p].size(); ++i) {
				data_pts.at<float>(c, 0) = pts[p][i].x;
				data_pts.at<float>(c, 1) = pts[p][i].y;
				--c;
			}
		}

		//Perform PCA analysis
		cv::PCA pca_analysis(data_pts, cv::Mat(), cv::PCA::DATA_AS_ROW, 2);

		//Store the center of the object
		ellipse.center = vec2f(pca_analysis.mean.at<float>(0, 0), pca_analysis.mean.at<float>(0, 1));

		//Store the eigenvalues and eigenvectors
		std::vector<vec2f> eigen_vecs(2);
		std::vector<float> eigen_val(2);
		for (int i = 0; i < 2; ++i)
		{
			eigen_vecs[i] = vec2f(pca_analysis.eigenvectors.at<float>(i, 0),
				pca_analysis.eigenvectors.at<float>(i, 1));
			//eigen_val[i] = pca_analysis.eigenvalues.at<float>(0, i); // Wrong values
		}
		assert(!dot(eigen_vecs[0], eigen_vecs[1]) && "Eigen vectors aren't orthogonal!");

		//find eigen values by projecting the points on eigenVectors
		findEigenValue(data_pts, eigen_vecs, eigen_val);

		//assign to R1
		for (int i = 0; i < 4; ++i) {
			const int indx = i % 2;
			R1(cs, i) = eigen_val[indx];
		}
		//assign to R2
		for (int i = 0; i < 4; ++i) {
			const int indx = (i + 1) % 2;
			R2(cs, i) = eigen_val[indx];
		}
		//assign to theta
		float angle = atan2(eigen_vecs[0].y, eigen_vecs[0].x);
		angle = angle < 0 ? angle + 2.f*pi : angle;

		theta(cs, 0) = angle;

		for (int i = 1; i < 4; ++i) {
			angle += pi / 2.f;
			angle = angle > 2.f*pi ? angle - 2.f*pi : angle;
			theta(cs, i) = angle;
		}

		//assign to isValid
		if (std::abs(eigen_val[0] - eigen_val[1]) < std::max(eigen_val[0], eigen_val[1]) / 2.f) {
			//std::cout << cs << " inValid is false. \n";
			isValid[cs] = false;
		}
	}
}

void CrossSection::greedyOpt(const Eigen::MatrixXf &R1, const Eigen::MatrixXf &R2, const Eigen::MatrixXf &theta, const std::vector<bool> &isValid, std::vector<Ellipse> &validEllipses) {
	const int plane_num = isValid.size();
	Ellipse ellps;
	ellps.longR = R1(0, 0);
	ellps.shortR = R2(0, 0);
	ellps.angle = theta(0, 0);
	validEllipses.push_back(ellps);
	float theta_old = ellps.angle;
	int opt_indx;
	for (int i = 1; i < plane_num; ++i) {
		if (!isValid[i]) {
			ellps.longR = validEllipses[i - 1].longR;
			ellps.shortR = validEllipses[i - 1].shortR;
			ellps.angle = validEllipses[i - 1].angle;
			validEllipses.push_back(ellps);
			continue;
		}
		//find the min theta with previous cross-section	
		float min_delta_theta = 2.f*pi;
		for (int j = 0; j < 4; ++j) {
			//argmin (theta_j - theta_old)
			float delta_theta = std::abs(theta(i, j) - theta_old);
			delta_theta = std::min(delta_theta, 2 * pi - delta_theta);

			if (delta_theta < min_delta_theta) {
				min_delta_theta = delta_theta;
				opt_indx = j;
			}
		}
		ellps.longR = R1(i, opt_indx);
		ellps.shortR = R2(i, opt_indx);
		ellps.angle = theta(i, opt_indx);
		theta_old = ellps.angle;
		validEllipses.push_back(ellps);
	}
}

void CrossSection::extractCompressParam(const std::vector<yarnIntersect2D> &allPlaneIntersect, std::vector<Ellipse> &ellipses) {

	//1. extract the ellipses 
	const int n = allPlaneIntersect.size();
	std::vector<Ellipse> allEllipses(n);
	std::vector<bool> isValid;
	fitEllipses(allPlaneIntersect, allEllipses, isValid);

	//2. precompute R1\R2\theta and isValid
	Eigen::MatrixXf R1(n, 4);
	Eigen::MatrixXf R2(n, 4);
	Eigen::MatrixXf theta(n, 4);
	preComputeEllipses(allEllipses, R1, R2, theta);
	
	//3. find optimal ellipse for valid ones
	std::vector<Ellipse> validEllipses;
	greedyOpt(R1, R2, theta, isValid, validEllipses);

	//4. interpolate for not valid ones
	std::vector<Ellipse> interEllipses;
	interpolateEllipses(validEllipses, isValid, interEllipses);

	//5. regularize ellipses
	//6. update R1, R2 and theta

	ellipses = interEllipses;


	std::cout << "Compression parameters for each cross-sections are written to the file! \n";
}

void CrossSection::retreiveSol(const Eigen::MatrixXf &R1, const Eigen::MatrixXf &R2, const Eigen::MatrixXf &theta, const Eigen::MatrixXf &totalCost, 
	const Eigen::MatrixXf &preConfig, const std::vector<bool> &isValid, std::vector<Ellipse> &validEllipses) {
	const int n = isValid.size();
	validEllipses.resize(n);
	std::vector<int> solutions(n);
	std::vector<float> cost_n{ totalCost(n-1,0), totalCost(n-1,1), totalCost(n-1,2), totalCost(n-1,3) };
	int opt_config = std::min_element(cost_n.begin(), cost_n.end()) - cost_n.begin();
	Ellipse ell;
	int i = n - 1;
	int ip = 0;
	for (int i = n -1; i > 0; --i )
	{
		if (R1(i, opt_config) != R1(i, opt_config) || R2(i, opt_config) != R2(i, opt_config) || theta(i, opt_config) != theta(i, opt_config)) {
			//std::cout << "\n" << i << "********************** " << opt_config << " " << R1(i, opt_config) << " " << R2(i, opt_config) << " " << theta(i, opt_config) << "\n";

			solutions[i] = 0;
			ell.longR = 1.f;
			ell.shortR = 0.f;
			ell.angle = 0.f;
			validEllipses[i] = ell;
		}
		else {
			solutions[i] = opt_config;
			ell.longR = R1(i, opt_config);
			ell.shortR = R2(i, opt_config);
			ell.angle = theta(i, opt_config);
			validEllipses[i] = ell;
		}
		
		opt_config = preConfig(i, opt_config);
		//if (isValid[i]) {
		//	std::cout << " " << i;
		//	ip = i - 1;
		//	while (!isValid[ip])
		//		ip = ip - 1;
		//	i = ip;
		//	std::cout << " **** " << i;
		//}
		//else
		//	i = i - 1;
	}
	ell.longR = R1(0, opt_config);
	ell.shortR = R2(0, opt_config);
	ell.angle = theta(0, opt_config);
	validEllipses[0] = ell;
}
void CrossSection::optimizeEllipses() {
	/****************************************************/

	////1. precompute R1\R2\theta and isValid
	//const int n = allPlaneIntersect.size();
	//Eigen::MatrixXf R1(n, 4);
	//Eigen::MatrixXf R2(n, 4);
	//Eigen::MatrixXf theta(n, 4);
	//std::vector<bool> isValid;
	//preComputeEllipses(allPlaneIntersect, R1, R2, theta, isValid);

	////2. precompute cost function 
	//std::vector<Eigen::Matrix4f> cost;
	//costFunction(R1, R2, theta, isValid, cost);

	////3. dynamic programming
	//Eigen::MatrixXf totalCost(n, 4); 
	//Eigen::MatrixXf preConfig(n, 4);
	//dynamicProgramming(isValid, cost, totalCost, preConfig);

	////4.retreive solution for valid cross-sections
	//std::vector<Ellipse> validEllipses(n);
	//std::vector<int> solutions(n);
	//std::vector<float> cost_n{ totalCost(n-1,0), totalCost(n-1,1), totalCost(n-1,2), totalCost(n-1,3) };
	//int opt_config = std::min_element(cost_n.begin(), cost_n.end()) - cost_n.begin();
	//Ellipse ell;
	//int i = n - 1;
	//int ip = 0;
	//while(i){

	//	solutions[i] = opt_config;
	//	ell.longR = R1(i, opt_config);
	//	ell.shortR = R2(i, opt_config);
	//	ell.angle = theta(i, opt_config);
	//	validEllipses[i] = ell;

	//	opt_config = preConfig(i, opt_config);
	//	if (isValid[i]) {
	//		ip = i - 1;
	//		while (!isValid[ip])
	//			ip = ip - 1;
	//		i = ip;
	//	}
	//	else
	//		i = i - 1;
	//}
	//ell.longR = R1(0, opt_config);
	//ell.shortR = R2(0, opt_config);
	//ell.angle = theta(0, opt_config);
	//validEllipses[0] = ell;

	////5. interpolate for not valid ones
	//interpolateEllipses(validEllipses, isValid, ellipses);
}

void CrossSection::parameterizeEllipses(const std::vector<Ellipse> &ellipses, std::vector<Ellipse> &simple_ellipses) {

	const int ignorPlanes = 0.15 * ellipses.size();
	float ell_lng_avg = 0.f, ell_shrt_avg = 0.f, ell_angle_avg = 0.f;
	float ell_shrt_min = std::numeric_limits<float>::max();
	float ell_shrt_max = std::numeric_limits<float>::min();
	for (int i = ignorPlanes; i < ellipses.size() - ignorPlanes; ++i)
	{
		ell_lng_avg += ellipses[i].longR;
		ell_shrt_avg += ellipses[i].shortR;
		ell_angle_avg += std::abs(ellipses[i].angle);
		if (ellipses[i].shortR < ell_shrt_min)
			ell_shrt_min = ellipses[i].shortR;
		if (ellipses[i].shortR > ell_shrt_max)
			ell_shrt_max = ellipses[i].shortR;
	}
	ell_lng_avg /= (ellipses.size() - 2 * ignorPlanes);
	ell_shrt_avg /= (ellipses.size() - 2 * ignorPlanes);
	ell_angle_avg /= (ellipses.size() - 2 * ignorPlanes);

	float ell_lng_avg_clusterL = 0.f;
	int ell_lng_avg_clusterL_num = 0;
	float ell_shrt_avg_clusterS = 0.f;
	int ell_shrt_avg_clusterS_num = 0;
	for (int i = ignorPlanes; i < ellipses.size() - ignorPlanes; ++i)
	{
		if (ellipses[i].longR > ell_lng_avg) {
			ell_lng_avg_clusterL += ellipses[i].longR;
			ell_lng_avg_clusterL_num++;
		}
		if (ellipses[i].shortR < ell_shrt_avg) {
			ell_shrt_avg_clusterS += ellipses[i].shortR;
			ell_shrt_avg_clusterS_num++;
		}
	}
	ell_lng_avg_clusterL /= static_cast<float>(ell_lng_avg_clusterL_num);
	ell_shrt_avg_clusterS /= static_cast<float>(ell_shrt_avg_clusterS_num);


	unsigned int numberOfPoints = ellipses.size();
	Point2DVector points;
	for (int i = ignorPlanes; i < ellipses.size() - ignorPlanes; ++i) {
		Eigen::Vector2d point;
		point(0) = i;
		point(1) = ellipses[i].shortR;
		points.push_back(point);
	}
	Eigen::VectorXd x(3);
	//x.fill(1.f);
	x(0) = (ell_shrt_max - ell_shrt_min) / 2.f;
	x(1) = 0.125; // 2.0*pi * 1.0 / 40.0;
	x(2) = (ell_shrt_max - ell_shrt_min) / 2.f;

	MyFunctorNumericalDiff functor;
	functor.Points = points;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);
	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(x);



	for (int i = 0; i < ellipses.size(); ++i) {
		Ellipse ell;

		ell.longR = ellipses[i].longR;
		ell.shortR = ellipses[i].shortR;
		ell.angle = ellipses[i].angle;

		//ell.longR = ell_lng_avg;
		//ell.longR = ell_lng_avg_clusterL;

		//ell.shortR = ell_shrt_avg;
		//ell.shortR = ell_shrt_avg_clusterS;
		//ell.shortR = x(0) * sin(x(1)*i) + x(2);

		//ell.angle = ell_angle_avg;

		simple_ellipses.push_back(ell);
	}
}

void CrossSection::planeIts2world(Plane &plane, vec2f &plane_point, vec3f &world_point) {
	//Finding 3D coord of a point of the plane having its 2D: inverse of project2plane()
	vec3f n = plane.n;
	vec3f e1 = plane.e1;
	vec3f e2 = plane.e2;

	vec3f local(plane_point.x, plane_point.y, 0.f); //as the point exists on the plane

	/*world = [e1 e2 n] * local
	local = [e1; e2; n] * world*/


	// ****** e1 maps to y and e2 maps to x *********
	//*************************************************
	world_point.x = dot(vec3f(e1.x, e2.x, n.x), local) + plane.point.x;
	world_point.y = dot(vec3f(e1.y, e2.y, n.y), local) + plane.point.y;
	world_point.z = dot(vec3f(e1.z, e2.z, n.z), local) + plane.point.z;
	//note that order of the axis matches with the project2plane()
}

#if 0
void CrossSection::extractNormals(std::vector<Ellipse> &ellipses, std::vector<vec3f> &normals, const char* normsFile) {

	FILE *fout;
	if (fopen_s(&fout, normsFile, "wt") == 0) {
		fprintf_s(fout, "%d\n", m_planesList.size());
		for (int i = 0; i < m_planesList.size(); ++i) {
			//rotate plane.e1=[1,0] by theta in the plane to obtain ellipse-short-axis
			vec2f p2D(cos(ellipses[i].angle), sin(ellipses[i].angle));
			//now project it to world coord
			vec3f end3D;
			planeIts2world(m_planesList[i], p2D, end3D);
			vec3f start3D = m_planesList[i].point;
			//note that ellipse center is same as plane.point
			vec3f normal = end3D - start3D;
			// vec3f normal = m_planesList[i].e1; clearly not the case
			// vec3f normal(1.f, 0.f, 0.f); clearly shaky
			assert(nv::length(normal) - 1.f < HERMITE_EPS   && "Normal vector is not normalized!");
			normals.push_back(normal);
			fprintf_s(fout, "%.6f %.6f %.6f \n", normal.x, normal.y, normal.z);
		}
	}
	fclose(fout);
}
#endif


#if 0
/* Given a yarn dataStructure, transform it to a vector of cross-sections (for debug use)*/
void CrossSection::yarn2crossSections(std::vector<yarnIntersect2D> &itsLists) {
	//first initialize the vectors
	itsLists.resize(m_planesList.size());
	for (int i = 0; i < itsLists.size(); ++i)
		itsLists[i].resize(m_yarn.plys.size());

	//copy the yarn into new dataStructure
	for (int p = 0; p < m_yarn.plys.size(); ++p) {
		plyItersect plyIts;
		for (int f = 0; f < m_yarn.plys[p].fibers.size(); ++f) {
			for (int v = 0; v < m_yarn.plys[p].fibers[f].vertices.size(); ++v) {
				itsLists[v][p].push_back(vec2f(m_yarn.plys[p].fibers[f].vertices[v].x,
					m_yarn.plys[p].fibers[f].vertices[v].y));
			}
		}
	}
}
#endif

void CrossSection::deCompressYarn(const std::vector<yarnIntersect2D> &planeIts, const float yarn_radius, std::vector<Ellipse> &ellipses, std::vector<yarnIntersect2D> &deCompressPlaneIts) {
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	deCompressPlaneIts.resize(planeIts.size());
	for (int i = 0; i < planeIts.size(); ++i) // number of planes
	{
		deCompressPlaneIts[i].resize(planeIts[i].size());
		for (int p = 0; p < planeIts[i].size(); ++p) // number of plys
		{
			deCompressPlaneIts[i][p].resize(planeIts[i][p].size());
			for (int v = 0; v < planeIts[i][p].size(); ++v) // number of intersections
			{
				vec2f its = planeIts[i][p][v];

				const float ellipse_long = ellipses[i].longR;
				const float ellipse_short = ellipses[i].shortR;
				const float ellipse_theta = ellipses[i].angle;

				//obtain ellipse axis
				const vec2f ellipse_axis_long = vec2f(cos(ellipse_theta), sin(ellipse_theta));
				const vec2f ellipse_axis_short = vec2f(-sin(ellipse_theta), cos(ellipse_theta));
				//project points on ellipse axis
				float _p_long = nv::dot(its, ellipse_axis_long);
				float _p_short = nv::dot(its, ellipse_axis_short);
				//apply the decompression
				_p_long *= yarn_radius / ellipse_long;
				_p_short *= yarn_radius / ellipse_short;
				vec2f _p_ellipse(_p_long, _p_short);
				//go back to e1-e2 space from ellipse space
				float _p_x = nv::dot(vec2f(ellipse_axis_long.x, ellipse_axis_short.x), _p_ellipse);
				float _p_y = nv::dot(vec2f(ellipse_axis_long.y, ellipse_axis_short.y), _p_ellipse);

				deCompressPlaneIts[i][p][v] = vec2f(_p_x, _p_y);
			}
		}
	}
}

void CrossSection::transferLocal2XY(const std::vector<yarnIntersect2D> &e1e2_Its, std::vector<yarnIntersect2D> &xy_Its) {
	//initialized the vector with a copy
	xy_Its = e1e2_Its;
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 0; i < e1e2_Its.size(); ++i) { // all planes
		vec3f n = m_planesList[i].n;
		vec3f e1 = m_planesList[i].e1;
		vec3f e2 = m_planesList[i].e2;
		int ply_num = e1e2_Its[i].size();
		for (int p = 0; p < ply_num; ++p) { // all plys
			vec2f plyCntr(0.f);
			for (int v = 0; v < e1e2_Its[i][p].size(); ++v) { // ply intersections
				vec2f e1e2_p = e1e2_Its[i][p][v];

				/* transform plyCenters in e1-e2 coord to xy plane */
				/****************************************************/
				// project e2 to ex, and e1 to ey with no translation
				xy_Its[i][p][v].x = dot(vec2f(e1.x, e2.x), e1e2_p);
				xy_Its[i][p][v].y = dot(vec2f(e1.y, e2.y), e1e2_p);
			}
		}
	}
}

void CrossSection::extractPlyTwist(const std::vector<yarnIntersect2D> &decompressPlaneIts, const char *plyCenterFile) {
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	std::ofstream fout(plyCenterFile);
	for (int i = 0; i < decompressPlaneIts.size(); ++i) { // all planes	
		int ply_num = decompressPlaneIts[i].size();
		for (int p = 0; p < ply_num; ++p) { // all plys
			vec2f plyCntr(0.f);
			for (int v = 0; v < decompressPlaneIts[i][p].size(); ++v) { // ply intersections
				plyCntr += decompressPlaneIts[i][p][v];
			}
			plyCntr /= static_cast<float>(decompressPlaneIts[i][p].size());

			fout << plyCntr.x << " " << plyCntr.y << '\n';
		}
		fout << '\n';
	}
	fout.close();
}

void CrossSection::parameterizePlyCenter(const char *plyCenterFile, const char *ParameterizePlyCntrFile) {
	std::vector<float> allR;
	std::vector<float> allTheta;
	std::ifstream fin(plyCenterFile);
	assert(!fin.fail());
	for (int i = 0; i < m_planesList.size(); ++i) {

		vec2f plyCntr1(0.f), plyCntr2(0.f);
		fin >> plyCntr1.x >> plyCntr1.y;
		fin >> plyCntr2.x >> plyCntr2.y;

		//find R between yarnCenter and plyCenter[0]
		float R = length(plyCntr1);
		float theta = atan2(plyCntr1.y, plyCntr1.x) + pi; //map between 0 to 2pi (NOTE that we should shift it back at the end ##)
		allR.push_back(R);
		allTheta.push_back(theta);
	}
	fin.close();


	/// regularization // *****************************
	std::vector<float> allR_regularized;
	std::vector<float> allTheta_regularized;
	allR_regularized.resize(allR.size());
	allTheta_regularized.resize(allTheta.size());
	//copy ellipse to simple_ellipse to fill up the two ends before and after sigma 
	for (int i = 0; i <  m_planesList.size(); ++i) {
		allR_regularized[i] = allR[i];
		allTheta_regularized[i] = allTheta[i];
	}

	const int sigma = 20;
	for (int i = sigma; i < m_planesList.size() - sigma; ++i) {
		float sum_R = allR[i];
		float sum_theta = allTheta[i]; //not valid for theta!!!

		for (int s = 0; s < sigma; ++s) {
			sum_R += allR[i - s];
			sum_R += allR[i + s];

			sum_theta += allTheta[i - s];
			sum_theta += allTheta[i + s];

		}
		float span = 2.f * sigma + 1.f;
		allR_regularized[i] = sum_R / span;
		allTheta_regularized[i] = sum_theta / span;
	}

	//////////// *****************************

	bool clockwise = false;
	const int ignorPlanes = 0.15 *  m_planesList.size();
	std::vector<int> periods;
	float jump_prev = ignorPlanes; //initialize with starting planeID
								   // to find when jumps are happened
	for (int i = ignorPlanes; i < m_planesList.size() - ignorPlanes; ++i) {
		float deltaTheta = std::abs(allTheta[i] - allTheta[i + 1]);
		if (deltaTheta > 2 * pi - deltaTheta) {
			periods.push_back(i - jump_prev);
			jump_prev = i;
			if (allTheta[i] - allTheta[i + 1] > 0)
				clockwise = 0; //theta gets larger after the jump for ccw
			else
				clockwise = 1;
		}
	}
	int period_avg = 0;
	for (int p = 1; p < periods.size(); p++)  //we ignor first period because they might not be a complete cycle
		period_avg += periods[p];
	period_avg /= (periods.size() - 1); //because the first cycle is discarded

	float R_avg = 0.f;
	float R_min = std::numeric_limits<float>::max();
	float R_max = std::numeric_limits<float>::min();
	for (int i = ignorPlanes; i < m_planesList.size() - ignorPlanes; ++i) {
		R_avg += allR[i];
		if (allR[i] > R_max)
			R_max = allR[i];
		if (allR[i] < R_min)
			R_min = allR[i];
	}

	R_avg /= static_cast<float>(m_planesList.size() - 2 * ignorPlanes);
	//now split the point as "below the avg" and "above" and get the avg of the first cluster
	//Because we are using fix ellipse params and ply-centrs get above the average for ellipse short
	float R_avg_clusterS = 0.f;
	int avg_cluster_numS = 0;
	float R_avg_clusterL = 0.f;
	int avg_cluster_numL = 0;
	for (int i = ignorPlanes; i < m_planesList.size() - ignorPlanes; ++i) {
		if (allR[i] < R_avg) {
			R_avg_clusterS += allR[i];
			avg_cluster_numS++;
		}
		else {
			R_avg_clusterL += allR[i];
			avg_cluster_numL++;
		}
	}
	R_avg_clusterS /= static_cast<float>(avg_cluster_numS);
	R_avg_clusterL /= static_cast<float>(avg_cluster_numL);

	std::ofstream fout(ParameterizePlyCntrFile);
	assert(!fout.fail());
	fout << m_planesList.size() << '\n';
	float theta;
	for (int i = 0; i < allR.size(); ++i) {
		if (clockwise)
			theta = 2 * pi - (i % period_avg) * 2.f * pi / period_avg;
		else
			theta = (i % period_avg) * 2.f * pi / period_avg;
		//TODO: subtract pi because we had shifted all by pi in ##

		//const float mag = (R_avg - R_avg_clusterS) / 2.f;
		//const float R_ = mag * sin(2.f * (theta - pi / 4.f - 1.5697)) + mag + R_avg_clusterS;

		//fout << allR_regularized[i] << " " << theta << '\n';
		fout << allR[i] << " " << allTheta[i] << '\n';
		//fout << R_avg_clusterS << " " << theta << '\n';
		//fout << R_avg << " " << theta << '\n'; 
	}
	fout.close();
}

void CrossSection::extractFiberVectors(const std::vector<yarnIntersect2D> &decompressPlaneIts, std::vector<plyItersect2D> &fiberCntrVector) {
	fiberCntrVector.resize(decompressPlaneIts.size());
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 0; i < decompressPlaneIts.size(); ++i) { // all planes	
		int p = 0; //since fiber-centers in only ply[0] are stable in e1-e2 space 
		vec2f plyCntr(0.f);
		for (int v = 0; v < decompressPlaneIts[i][p].size(); ++v) { // ply intersections
			plyCntr += decompressPlaneIts[i][p][v];
		}
		plyCntr /= static_cast<float>(decompressPlaneIts[i][p].size());
		//Find fiber[0]-plyCenter for each plane
		for (int v = 0; v < decompressPlaneIts[i][p].size(); ++v) {
			vec2f fiberVec = decompressPlaneIts[i][p][v] - plyCntr;
			fiberCntrVector[i].push_back(fiberVec);
		}
	}
}
void CrossSection::fiberTwisting(const std::vector<yarnIntersect2D> &decompressPlaneIts, std::vector<float> &fiber_theta, const char *fiberTwistFile) {

	std::vector<plyItersect2D> fiberCntrVector;
	extractFiberVectors(decompressPlaneIts, fiberCntrVector);

	//push theta = 0 for the first plane
	fiber_theta.push_back(0.f);

	std::ofstream fout(fiberTwistFile);
	fout << fiberCntrVector.size() << '\n';
	fout << fiber_theta[0] << '\n';
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 0; i < fiberCntrVector.size() - 1; ++i) { // all planes (-1 because it checks the angle between two consecutive points)
		float theta_avg = 0;
		const int fiber_num = std::min(fiberCntrVector[i].size(), fiberCntrVector[i + 1].size()); //At the endpoints, next plane might not have as many as the current plane
		for (int v = 0; v < fiber_num; ++v) { // fiber-centers for each plane
			//compute the angle between two vec2f
			float cos = dot(fiberCntrVector[i][v], fiberCntrVector[i + 1][v]) / (length(fiberCntrVector[i + 1][v]) * length(fiberCntrVector[i][v]));
			if (1.f - cos > 1e-7) // skip the nan that are generated duo to small acos()
				theta_avg += acos(cos);
		}

		theta_avg /= fiber_num;
		fiber_theta.push_back(theta_avg);
		fout << fiber_theta[i + 1] << '\n';
	}
	fout.close();
}


void CrossSection::regularizeEllipses(const std::vector<Ellipse> &ellipses, const std::vector<float> &theta_R, std::vector<Ellipse> &simple_ellipses, std::vector<float> &simple_theta_R, const int sigma) {
	
	simple_theta_R.resize(ellipses.size());
	simple_ellipses.resize(ellipses.size());
	//copy ellipse to simple_ellipse to fill up the two ends before and after sigma 
	for (int i = 0; i < ellipses.size() ; ++i) {
		simple_ellipses[i].longR = ellipses[i].longR;
		simple_ellipses[i].shortR = ellipses[i].shortR;
		simple_ellipses[i].angle = ellipses[i].angle;
		simple_theta_R[i] = theta_R[i];
	}

	//gaussian smoothing
	for (int i = sigma; i < ellipses.size() - sigma; ++i) {
		float sum_lng = 0.f;
		float sum_shrt = 0.f;
		float sum_angle = 0.f;
		float sum_theta_R = 0.f;

		//find if in the sigma terithory, there is any fluctuation
		bool isFluctuate = false;
		for (int s = 0; s < sigma; ++s) {
			float delta = std::abs(ellipses[i - s].angle - ellipses[i].angle);
			if (delta > pi) //since we don't normally have this rotation, this is absolutely the case of fluctuation at 2pi-0
			{
				isFluctuate = true;
				break;
			}
		}
		bool isFluctuate_R = false;
		for (int s = 0; s < sigma; ++s) {
			float delta = std::abs(theta_R[i - s] - theta_R[i]);
			if (delta > pi) //since we don't normally have this rotation, this is absolutely the case of fluctuation at 2pi-0
			{
				isFluctuate_R = true;
				break;
			}
		}

		for (int s = 0; s < sigma; ++s) {
			sum_lng += ellipses[i - s].longR;
			sum_lng += ellipses[i + s].longR;

			sum_shrt += ellipses[i - s].shortR;
			sum_shrt += ellipses[i + s].shortR;

			float angle_s = ellipses[i - s].angle;
			float angle__s = ellipses[i + s].angle;
			//handle the case of 2*pi and zero
			if (isFluctuate) { //negate angle[i-s] if it's on the wrong side
				angle_s = ellipses[i - s].angle > pi ? ellipses[i - s].angle - 2.f * pi : ellipses[i - s].angle;
				angle__s = ellipses[i + s].angle > pi ? ellipses[i + s].angle - 2.f * pi : ellipses[i + s].angle;
			}			
			sum_angle += angle_s;
			sum_angle += angle__s;

			/* theta_R */
			float angle_s_R = theta_R[i-s];
			float angle__s_R = theta_R[i + s];
			//handle the case of 2*pi and zero
			if (isFluctuate_R) { //negate angle[i-s] if it's on the wrong side
				angle_s_R = theta_R[i - s] > pi ? theta_R[i - s] - 2.f * pi : theta_R[i - s];
				angle__s_R = theta_R[i + s] > pi ? theta_R[i + s] - 2.f * pi : theta_R[i + s];
			}
			sum_theta_R += angle_s_R;
			sum_theta_R += angle__s_R;
		}

		sum_lng += ellipses[i].longR;
		sum_shrt += ellipses[i].shortR;
		sum_angle += ( isFluctuate && ellipses[i].angle > pi ) ? ellipses[i].angle - 2.f * pi : ellipses[i].angle;
		sum_theta_R += (isFluctuate_R && theta_R[i] > pi) ? theta_R[i] - 2.f * pi : theta_R[i];

		float span = 2.f * sigma + 1.f;
		simple_ellipses[i].longR = sum_lng / span;
		simple_ellipses[i].shortR = sum_shrt / span;
		float avg_angle = sum_angle / span;
		simple_ellipses[i].angle = avg_angle < 0 ? 2.f*pi + avg_angle : avg_angle;

		float avg_theta_R = sum_theta_R / span;
		simple_theta_R[i] = avg_theta_R < 0 ? 2.f*pi + avg_theta_R : avg_theta_R;

		//simple_ellipses[i].longR = ellipses[i].longR;
		//simple_ellipses[i].shortR = ellipses[i].shortR;
		//simple_ellipses[i].angle = ellipses[i].angle;

	}
	//debug only:
	std::ofstream fout("compress_regularized.txt");
	fout << ellipses.size() << '\n';
	for (int i = 0; i < ellipses.size(); ++i) {
		fout << simple_ellipses[i].longR << " " << simple_ellipses[i].shortR << " " << simple_ellipses[i].angle << '\n';
	}
	fout.close();
}


void CrossSection::regularizeEllipses(const std::vector<Ellipse> &ellipses, std::vector<Ellipse> &simple_ellipses, const int sigma) {

	simple_ellipses.resize(ellipses.size());
	//copy ellipse to simple_ellipse to fill up the two ends before and after sigma 
	for (int i = 0; i < ellipses.size(); ++i) {
		simple_ellipses[i].longR = ellipses[i].longR;
		simple_ellipses[i].shortR = ellipses[i].shortR;
		simple_ellipses[i].angle = ellipses[i].angle;
	}

	//gaussian smoothing
	for (int i = sigma; i < ellipses.size() - sigma; ++i) {
		float sum_lng = 0.f;
		float sum_shrt = 0.f;
		float sum_angle = 0.f;

		//find if in the sigma terithory, there is any fluctuation
		bool isFluctuate = false;
		for (int s = 0; s < sigma; ++s) {
			float delta = std::abs(ellipses[i - s].angle - ellipses[i].angle);
			if (delta > pi) //since we don't normally have this rotation, this is absolutely the case of fluctuation at 2pi-0
			{
				isFluctuate = true;
				break;
			}
		}

		for (int s = 0; s < sigma; ++s) {
			sum_lng += ellipses[i - s].longR;
			sum_lng += ellipses[i + s].longR;

			sum_shrt += ellipses[i - s].shortR;
			sum_shrt += ellipses[i + s].shortR;

			float angle_s = ellipses[i - s].angle;
			float angle__s = ellipses[i + s].angle;
			//handle the case of 2*pi and zero
			if (isFluctuate) { //negate angle[i-s] if it's on the wrong side
				angle_s = ellipses[i - s].angle > pi ? ellipses[i - s].angle - 2.f * pi : ellipses[i - s].angle;
				angle__s = ellipses[i + s].angle > pi ? ellipses[i + s].angle - 2.f * pi : ellipses[i + s].angle;
			}
			sum_angle += angle_s;
			sum_angle += angle__s;
		}

		sum_lng += ellipses[i].longR;
		sum_shrt += ellipses[i].shortR;
		sum_angle += (isFluctuate && ellipses[i].angle > pi) ? ellipses[i].angle - 2.f * pi : ellipses[i].angle;

		float span = 2.f * sigma + 1.f;
		simple_ellipses[i].longR = sum_lng / span;
		simple_ellipses[i].shortR = sum_shrt / span;
		float avg_angle = sum_angle / span;
		simple_ellipses[i].angle = avg_angle < 0 ? 2.f*pi + avg_angle : avg_angle;

		//simple_ellipses[i].longR = ellipses[i].longR;
		//simple_ellipses[i].shortR = ellipses[i].shortR;
		//simple_ellipses[i].angle = ellipses[i].angle;

	}
	//debug only:
	std::ofstream fout("compress_regularized.txt");
	fout << ellipses.size() << '\n';
	for (int i = 0; i < ellipses.size(); ++i) {
		fout << simple_ellipses[i].longR << " " << simple_ellipses[i].shortR << " " << simple_ellipses[i].angle << '\n';
	}
	fout.close();
}