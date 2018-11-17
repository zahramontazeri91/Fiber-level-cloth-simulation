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
	PlanesIntersections2D(itsLists, allPlaneIntersect);
}
void CrossSection::init_norm(const char* yarnfile, const int ply_num, const char* curvefile, const char* normfile,
	const int seg_subdiv, const int num_planes, std::vector<yarnIntersect2D> &allPlaneIntersect) {
	m_yarn.build(yarnfile, ply_num);
	m_curve.init_norm(curvefile, normfile, seg_subdiv);
	std::vector<yarnIntersect> itsLists;
	buildPlanes(num_planes, itsLists);
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

void CrossSection::shapeMatch_A(const Eigen::MatrixXf &pnt_trans, const Eigen::MatrixXf &pnt_ref, Eigen::Matrix2f &A, Eigen::Matrix2f &R) {

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

	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
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
	R = U*V.transpose();
	theta_R = atan2(R(1, 0), R(0, 0));

	//Return only S for generating training because global-rotation (R) may distract NN but store it for phase-matching purpose
	A = S;

}

void convertYarnIntersect2Mat(const yarnIntersect2D &pnts_yarn, Eigen::MatrixXf &pnts_mat) {
	//Find the total number of points for all plys
	int sz = 0;
	for (int p = 0; p < pnts_yarn.size(); ++p)
		sz += pnts_yarn[p].size();

	pnts_mat = Eigen::MatrixXf(2, sz);

	int c = 0;
	for (int p = 0; p < pnts_yarn.size(); ++p) {
		for (int i = 0; i < pnts_yarn[p].size(); ++i) {
			pnts_mat(0, c) = pnts_yarn[p][i].x;
			pnts_mat(1, c) = pnts_yarn[p][i].y;
			++c;
		}
	}
}

void writeFirstPnts(const std::vector<yarnIntersect2D> &all_pnts, const char* pnts_file) {
	const int cs_num = 1; //write only first cross-seciton
	std::ofstream fout(pnts_file);
	for (int c = 0; c < cs_num; c++) {
		const int plyNum = all_pnts[c].size();
		for (int p = 0; p < plyNum; p++) {
			const int fiberNum_ply = all_pnts[c][p].size();
			for (int f = 0; f < fiberNum_ply; f++) {
				fout << all_pnts[c][p][f].x << " " << all_pnts[c][p][f].y << std::endl; //total points are fibernum for each ply
			}
		}
	}
	fout.close();
	std::cout << "Average points are written to " << pnts_file << std::endl;
}
/*
void writeAvgPnts(const std::vector<yarnIntersect2D> &all_pnts, const int fiber_num, const char* pnts_file) {
	const int cs_num = all_pnts.size();
	vec2f pnt(0.0);
	std::vector<vec2f> avg(fiber_num, pnt);
	for (int c = 0; c < cs_num; c++) {
		const int plyNum = all_pnts[c].size();
		for (int p = 0; p < plyNum; p++) {
			const int fiberNum_ply = all_pnts[c][p].size();
			for (int f = 0; f < fiberNum_ply; f++) {
				avg[p*fiberNum_ply + f] = avg[p*fiberNum_ply + f] + all_pnts[c][p][f]; //total points are fibernum for each ply
			}
		}
	}
	std::ofstream fout(pnts_file);
	for (int i = 0; i < fiber_num; i++) {
		avg[i] = avg[i] / cs_num;
		fout << avg[i].x << " " << avg[i].y << std::endl;
	}
	fout.close();
	std::cout << "Average points are written to " << pnts_file << std::endl;
}


void writePnts(const std::vector<yarnIntersect2D> &all_pnts, const std::vector<Eigen::Matrix2f> &all_R, 
	const char* pnts_file, const int isRotate, const int ws_ds, const float trimPercent, const int sampleRate) {

	assert(all_R.size() == all_pnts.size());
	const int vrtx = all_R.size();

	const int ignorPlanes = trimPercent * vrtx; // crop the first and last #% of the yarn
	const int ws_us = ws_ds * sampleRate; // first upsample the window-size
	const int window_num = ((vrtx - ws_us + 1) - 2 * ignorPlanes);

	std::vector<int> indices;
	for (int w = ignorPlanes; w < (vrtx - ws_us + 1) - ignorPlanes; w = w + sampleRate)
		indices.push_back(w);

	std::ofstream fout(pnts_file);
	for (int omp_i = 0; omp_i < static_cast<int>(indices.size()); ++omp_i) {
		const int w = indices[omp_i];

		//define a curve segment 
		const int start = w;
		const int end = w + (ws_us - 1);
		const int v_yarn = ceil((start + end) / 2.0);

		Eigen::MatrixXf pnts_mat;
		convertYarnIntersect2Mat(all_pnts[v_yarn], pnts_mat);
		for (int j = 0; j < pnts_mat.cols(); j++) {
			if (isRotate) {
				Eigen::MatrixXf mat_rot(2, 1);
				mat_rot = all_R[v_yarn] * pnts_mat.col(j);
				fout << mat_rot(0, 0) << " " << mat_rot(1, 0) << std::endl;
			}
			else
				fout << pnts_mat(0, j) << " " << pnts_mat(1, j) << std::endl;
		}
		fout << std::endl;

	}
	fout.close();
}
*/


void CrossSection::yarnShapeMatch_A(const yarnIntersect2D &pnts_trans, const yarnIntersect2D &pnts_ref, Eigen::Matrix2f &A, Eigen::Matrix2f &R) {

	Eigen::MatrixXf all_ref;
	Eigen::MatrixXf all_trans;
	convertYarnIntersect2Mat(pnts_trans, all_trans);
	convertYarnIntersect2Mat(pnts_ref, all_ref);

	const int sz_ref = all_ref.cols();
	const int sz_trans = all_trans.cols();
	if (sz_ref != sz_trans) std::cout << "size ref-yarn and simulated are: " << sz_ref << " " << sz_trans << std::endl;
	assert(sz_ref == sz_trans);

	shapeMatch_A(all_trans, all_ref, A, R);
}


void CrossSection::yarnShapeMatches_A(const std::vector<yarnIntersect2D> &pnts_trans, const std::vector<yarnIntersect2D> &pnts_ref,
	std::vector<Eigen::Matrix2f> &all_A, std::vector<Eigen::Matrix2f> &all_R) {

	if (pnts_trans.size() != pnts_ref.size()) std::cout << pnts_trans.size() << " " << pnts_ref.size() << std::endl;
	assert(pnts_trans.size() == pnts_ref.size());
	const int n = pnts_trans.size();
	all_A.resize(n);
	all_R.resize(n);

	for (int i = 0; i < n; ++i) {
		Eigen::Matrix2f A, R;
		if (pnts_trans[i][0].size() != pnts_ref[i][0].size() || pnts_trans[i][1].size() != pnts_ref[i][1].size()) { //for two plys!
			std::cout << "Mismatch in " << i << "-th cross-section. fibersim: " << pnts_trans[i][0].size() << " ref: " << pnts_ref[i][0].size() << std::endl;
			A(0, 0) = 0.0;
			A(0, 1) = 1.0;
			A(1, 0) = 0.0;
			A(1, 1) = 1.0;
			R = A;
		}
		else {
			yarnShapeMatch_A(pnts_trans[i], pnts_ref[i], A, R);
		}
		all_A[i] = A;
		all_R[i] = R;
	}
}


void CrossSection::allPlyCenters(std::vector<std::vector<vec3f>> &plyCenters, std::vector<yarnIntersect> &itsLists) {

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






