#include "CrossSection.h"
#include "curveFitting.h"
#include <Eigen/SVD>

void CrossSection::init(const char* yarnfile, const int ply_num, const char* curvefile, const char* normfile,
	const int seg_subdiv, const int num_planes, std::vector<yarnIntersect2D> &allPlaneIntersect) {
	m_yarn.build(yarnfile, ply_num);
	//TODO: pass this to the func
	m_curve.init(curvefile, normfile, seg_subdiv);
	std::vector<yarnIntersect> itsLists;
	buildPlanes(num_planes, itsLists);
	std::cout << "Finding the intersections with the cross-sectional plane... \n";
	PlanesIntersections2D(itsLists, allPlaneIntersect);
}
void CrossSection::init_norm(const char* yarnfile, const int ply_num, const char* curvefile, const char* normfile,
	const int seg_subdiv, const int num_planes, std::vector<yarnIntersect2D> &allPlaneIntersect) {
	m_yarn.build(yarnfile, ply_num);
	m_curve.init_norm(curvefile, normfile, seg_subdiv);
	std::vector<yarnIntersect> itsLists;
	buildPlanes(num_planes, itsLists);
	std::cout << "Finding the intersections with the cross-sectional plane... \n";
	PlanesIntersections2D(itsLists, allPlaneIntersect);
}


void CrossSection::buildPlanes(const int num_planes, std::vector<yarnIntersect> &itsLists, const char* twistFile, const int upsample) { //TODO: merge it with next func
	
	m_planesList.resize(num_planes);
	std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;
	//m_curve.assign(all_pts, all_tang, all_norm);
	m_curve.assign_twist(twistFile, all_pts, all_tang, all_norm, upsample);
	assert(all_pts.size() == num_planes);

	for (int i = 0; i < num_planes; ++i) {

		Eigen::Vector3d ez = all_tang[i];
		Eigen::Vector3d ey = all_norm[i];
		Eigen::Vector3d ex = ez.cross(ey);
		m_planesList[i].point = vec3f(all_pts[i][0], all_pts[i][1], all_pts[i][2]);


		m_planesList[i].n = vec3f(ez[0], ez[1], ez[2]);
		m_planesList[i].e1 = vec3f(ex[0], ex[1], ex[2]);
		m_planesList[i].e2 = vec3f(ey[0], ey[1], ey[2]);
	}
	std::vector<std::vector<vec3f>> plyCenters;
	allPlanesIntersections(itsLists);
	allPlyCenters(plyCenters, itsLists);

}



void CrossSection::buildPlanes(const int num_planes, std::vector<yarnIntersect> &itsLists) {
	const double curveLength = m_curve.totalLength();
	m_planesList.resize(num_planes);

	std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;
	m_curve.assign(all_pts, all_tang, all_norm);
	assert(all_pts.size() == num_planes);

	for (int i = 0; i < num_planes; ++i) {

		Eigen::Vector3d ez = all_tang[i];
		Eigen::Vector3d ey = all_norm[i];
		Eigen::Vector3d ex = ez.cross(ey);
		m_planesList[i].point = vec3f(all_pts[i][0], all_pts[i][1], all_pts[i][2]);

		
		m_planesList[i].n = vec3f(ez[0], ez[1], ez[2]);
		m_planesList[i].e1 = vec3f(ex[0], ex[1], ex[2]);
		m_planesList[i].e2 = vec3f(ey[0], ey[1], ey[2]);

	}
	std::vector<std::vector<vec3f>> plyCenters;
	allPlanesIntersections(itsLists);

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

			float dist = length(its - plane.point);
			//assert(length(its - plane.point) > 1e-6 && "intersection exactly at the plane.point"); //TODO: handle this corner case later 
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
	bool isIntrsct = false;

	for (int i = 0; i < m_planesList.size(); ++i) {
		yarnIntersect itsList;
		if (yarnPlaneIntersection(m_planesList[i], itsList)) {
			isIntrsct = true;
			itsLists.push_back(itsList);
		}
		else {
			// if plane doesn't intersect the yarn, just add some values and will trim it afterward based on trimPercent
			vec3f tmp(0.f);
			plyIntersect tmp2;
			tmp2.push_back(tmp);
			itsList.push_back(tmp2);
			itsLists.push_back(itsList);
		}
	}
	if (m_planesList.size() != itsLists.size())
		std::cout << itsLists.size() << " out of " << m_planesList.size() << " many planes had intersections! \n";
	assert(m_planesList.size() == itsLists.size() && "not all the planes have intersections!");

	if (isIntrsct)
		return true;
	return false;
}

void CrossSection::shapeMatch_A(const Eigen::MatrixXf &pnt_trans, const Eigen::MatrixXf &pnt_ref, Eigen::Matrix2f &A, std::ofstream &phase_fout) {

	Matrix_S mat_S; float theta_R;

	assert(pnt_trans.cols() == pnt_ref.cols());
	const int n = pnt_trans.cols();
	Eigen::Matrix2f Apq = Eigen::Matrix2f::Zero();
	Eigen::Matrix2f Aqq_1 = Eigen::Matrix2f::Zero();
	Eigen::Matrix2f Aqq = Eigen::Matrix2f::Zero();

	Eigen::MatrixXf Qtrans = pnt_ref.transpose();
	Eigen::MatrixXf cm_ref = Eigen::MatrixXf::Zero(2, 1);
	Eigen::MatrixXf cm_trans = Eigen::MatrixXf::Zero(2, 1);
	if (0) { //allow translation
		for (int i = 0; i < n; ++i) {
			cm_ref += pnt_ref.col(i);
			cm_trans += pnt_trans.col(i);
		}
		cm_ref /= n;
		cm_trans /= n;
	}


	for (int i = 0; i < n; ++i) {
		Apq += (pnt_trans.col(i) - cm_trans)  *  (pnt_ref.col(i) - cm_ref).transpose() / n;
		Aqq_1 += (pnt_ref.col(i) - cm_ref) *  (pnt_ref.col(i) - cm_ref).transpose() / n;
	}

	assert(Aqq_1.determinant() != 0);
	Aqq = Aqq_1.inverse();
	A = Apq*Aqq;

	assert(Aqq_1.determinant() != 0);
	Aqq = Aqq_1.inverse();

	//SVD decomposition
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(Apq*Aqq, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::Matrix2f U = svd.matrixU();
	Eigen::Matrix2f sigma = svd.singularValues().asDiagonal();
	Eigen::Matrix2f V = svd.matrixV();

	//SR decompose
	Eigen::Matrix2f S = U*sigma*U.transpose();
	//we don't care about R at this point:
	Eigen::Matrix2f R = U*V.transpose();
	theta_R = atan2(R(1, 0), R(0, 0));
	phase_fout << R(0, 0) << " " << R(0, 1) << " " << R(1, 0) << " " << R(1, 1) << std::endl;

	//Return only S for generating training because global-rotation (R) may distract NN
	A = S;

}

void CrossSection::yarnShapeMatch_A(const yarnIntersect2D &pnts_trans, const yarnIntersect2D &pnts_ref, Eigen::Matrix2f &A, std::ofstream &phase_fout) {
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
	shapeMatch_A(all_trans, all_ref, A, phase_fout);
}


void CrossSection::yarnShapeMatches_A(const std::vector<yarnIntersect2D> &pnts_trans, const std::vector<yarnIntersect2D> &pnts_ref,
	std::vector<Eigen::Matrix2f> &all_A, std::ofstream &phase_fout) {

	assert(pnts_trans.size() == pnts_trans.size());
	const int n = pnts_trans.size();
	all_A.resize(n);
	const int num_of_cores = omp_get_num_procs();



	for (int i = 0; i < n; ++i) {
		Eigen::Matrix2f A;
		if (pnts_trans[i][0].size() != pnts_ref[i][0].size() || pnts_trans[i][1].size() != pnts_ref[i][1].size() ) { //for two plys!
			std::cout << "Not equal number of points in simulated and procedural in " << i << "-th cross-section.\n";
			std::cout << "trans: " << pnts_trans[i][0].size() << " ref: " << pnts_ref[i][0].size() << std::endl;
			A(0, 0) = 0.0;
			A(0, 1) = 1.0;
			A(1, 0) = 0.0;
			A(1, 1) = 1.0;
			phase_fout << 1.f << " " << 0.f << " " << 0.f <<  " " << 1.0 << std::endl; //set global rotation to 0-degree for lost planes
			continue;
		}
		yarnShapeMatch_A(pnts_trans[i], pnts_ref[i], A, phase_fout);
		all_A[i] = A;
	}
}


void CrossSection::allPlyCenters(std::vector<std::vector<vec3f>> &plyCenters, std::vector<yarnIntersect> &itsLists) {

	const int plane_num = itsLists.size();
	plyCenters.resize(plane_num); //number of planes
	const int num_of_cores = omp_get_num_procs();

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

	const int ply_num = m_yarn.plys.size();
	const int num_of_cores = omp_get_num_procs();


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

void interpolateEllipses(const std::vector<Ellipse> &ellipses, const std::vector<bool> &isValid, std::vector<Ellipse> &ellipses_smooth) {
	ellipses_smooth.resize(ellipses.size());
	assert(ellipses.size() == isValid.size());
	assert(isValid[0]); //assumption: handle this later
	int i = 0;
	while (i < isValid.size()) {
		if (!isValid[i]) { // i-1 is valid 
			int j = i + 1;
			while (!isValid[j]) // j is valid
				j++;
			//interpolate values between i-1 and j 
			float delta_indx = (j)-(i - 1) - 1;
			float delta_longR = ellipses[j].longR - ellipses[i - 1].longR;
			float delta_shortR = ellipses[j].shortR - ellipses[i - 1].shortR;
			float delta_angle = ellipses[j].angle - ellipses[i - 1].angle;
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
			isValid[cs] = false;
		}
	}
}

/* Given a yarn dataStructure, transform it to a vector of cross-sections (for debug use)*/
void CrossSection::yarn2crossSections(std::vector<yarnIntersect2D> &itsLists) {
	//first initialize the vectors
	int vrtx_num = m_yarn.getStepNum();
	itsLists.resize(vrtx_num);
	for (int i = 0; i < itsLists.size(); ++i) 
		itsLists[i].resize(m_yarn.plys.size());

	

	//copy the yarn into new dataStructure
	for (int p = 0; p < m_yarn.plys.size(); ++p) {
		plyIntersect plyIts;
		for (int f = 0; f < m_yarn.plys[p].fibers.size(); ++f) {

			for (int v = 0; v < m_yarn.plys[p].fibers[f].vertices.size(); ++v) {
				itsLists[v][p].push_back(vec2f(m_yarn.plys[p].fibers[f].vertices[v].x, m_yarn.plys[p].fibers[f].vertices[v].y));
			}
		}
	}
}






