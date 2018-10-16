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






