#include "CrossSection.h"
#include "curveFitting.h"

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

void CrossSection::init(const char* yarnfile, const int ply_num, const char* curvefile, const int seg_subdiv, const int num_planes) {
	m_yarn.build(yarnfile, ply_num);
	//std::cout << m_yarn.plys[0].fibers[10].vertices[10].x << std::endl;
	m_curve.init(curvefile, seg_subdiv);
	buildPlanes(num_planes);
}

void CrossSection::buildPlanes(const int num_planes) {
	const double curveLen = m_curve.totalLength();
	const double crossSectionLen = curveLen / static_cast<double>(num_planes - 1); //place plane at the very ends as well
	const double crossSection_t = m_curve.arcLengthInvApprox(crossSectionLen);
	m_planesList.resize(num_planes);
	for (int i = 0; i < num_planes; ++i) {
		Eigen::Vector3d curve_p = m_curve.eval(i * crossSection_t);
		Eigen::Vector3d curve_t = m_curve.evalTangent(i * crossSection_t);
		/* need binormal and normal of the plane to project 3D points on the plane */
		Eigen::Vector3d curve_n = m_curve.evalNormal(i * crossSection_t); //short axis
		Eigen::Vector3d curve_b = curve_t.cross(curve_n); //long axis
		m_planesList[i].point = vec3f(curve_p[0], curve_p[1], curve_p[2]);
		m_planesList[i].n     = vec3f(curve_t[0], curve_t[1], curve_t[2]);
		//m_planesList[i].e1    = vec3f(curve_n[0], curve_n[1], curve_n[2]);
		//m_planesList[i].e2 = cross(m_planesList[i].n, m_planesList[i].e1);
		//assert(dot(m_planesList[i].n, m_planesList[i].e1) < EPS);
	}
	std::vector<std::vector<vec3f>> plyCenters;
	allPlyCenters(plyCenters);
	// Define inplane 2D coord using the direction from yarn-center to intersection of first ply-center 
	for (int i = 0; i < num_planes; ++i) {
		// plyCenters[i][0] - m_planesList[i].point is too small so its dot product with n doesn't show they are perpendicular
		// std::cout << dot(e1, m_planesList[i].n) << std::endl;
		m_planesList[i].e1 = nv::normalize(plyCenters[i][0] - m_planesList[i].point);
		m_planesList[i].e2 = cross(m_planesList[i].n, m_planesList[i].e1);
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
		if (t >= EPS && t <= t_end) {
			its = start + dir*t;

			assert(length(its - plane.point) > 1e-6 && "intersection at the plane.point"); //TODO: handle this corner case later 
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

	assert( nv::length(plane.n) - 1.f < EPS   && "normal vector is not normalized!" );

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
				if (linePlaneIntersection(start, end, plane, its) ) {				

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

void CrossSection::allPlyCenters(std::vector<std::vector<vec3f>> &plyCenters) {

	std::vector<yarnIntersect> itsLists;
	allPlanesIntersections(itsLists);

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

	assert(m_planesList.size() == itsLists.size() );
		//std::cout << itsLists.size() << " out of " << m_planesList.size()  << " many planes had intersections! \n";
	
	const int ply_num = m_yarn.plys.size();
	FILE *fout;
	if (fopen_s(&fout, filename, "wt") == 0) {

		fprintf_s(fout, "plane_num: %d \n", itsLists.size());
		fprintf_s(fout, "ply_num: %d \n", ply_num);
		fprintf_s(fout, "\n");

		for (int cs = 0; cs < itsLists.size() ; ++cs) { //for each plane //TODO: why -1?
			/* First write the yarn-center */
			//fprintf_s(fout, "center: %.4lf %.4lf %.4lf \n", m_planesList[cs].point[0], m_planesList[cs].point[1], m_planesList[cs].point[2]);
			/* Then all the intersections for each ply */
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
	std::cout << "Intersections are written to the file successfully! \n\n";
}

void CrossSection::project2Plane (const vec3f& P3d, const Plane& plane, vec2f& P2d) {
	vec3f n = plane.n;
	vec3f e1 = plane.e1;
	vec3f e2 = plane.e2;
	
	assert(length(e1) - 1.f < 1e-6);
	assert(length(e2) - 1.f < 1e-6);
	assert(dot(e1, e2) < EPS);
	P2d.x = dot(e1, (P3d -  plane.point));
	P2d.y = dot(e2, (P3d - plane.point));
}

void CrossSection::PlanesIntersections2D(std::vector<yarnIntersect> &itsLists, std::vector<yarnIntersect2D> &allPlaneIntersect) {
	if (m_planesList.size() != itsLists.size())
		std::cout << itsLists.size() << " out of " << m_planesList.size() << " many planes had intersections! \n";

	const int ply_num = m_yarn.plys.size();

	for (int cs = 0; cs < itsLists.size(); ++cs) { //for each plane
		Plane plane;
		get_plane(cs, plane);
		/* First write the yarn-center */
		vec2f center2D;
		project2Plane(plane.point, plane, center2D);
		assert (center2D == vec2f(0.f) ) ;

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
	std::cout << "Intersections are written to the file successfully! \n\n";
}


void CrossSection::fitEllipse(const yarnIntersect2D &pts, Ellipse &ellipse)
{
	//Find the total number of points for all plys
	int sz = 0;
	for (int p = 0; p < pts.size(); ++p)
		sz += pts[p].size();
	cv::Mat data_pts (sz, 2, CV_32FC1, cv::Scalar::all(0));
	
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
	float max0 = std::numeric_limits<float>::min();
	float max1 = std::numeric_limits<float>::min();
	for (int i = 0; i < data_pts.rows; ++i) {
		float prj0 = eigen_vecs[0].x * data_pts.at<float>(i, 0) + eigen_vecs[0].y * data_pts.at<float>(i, 1);
		if ( std::abs(prj0) > max0)
			max0 = std::abs(prj0);
		float prj1 = eigen_vecs[1].x * data_pts.at<float>(i, 0) + eigen_vecs[1].y * data_pts.at<float>(i, 1);
		if (std::abs(prj1) > max1)
			max1 = std::abs(prj1);
	}
	// TODO: find the eigen values using eigen vectors
	eigen_val[0] = std::max(max0,max1);
	eigen_val[1] = std::min(max1,max0);

	vec2f p1 = ellipse.center + vec2f(eigen_vecs[0].x * eigen_val[0], eigen_vecs[0].y * eigen_val[0]);
	vec2f p2 = ellipse.center - vec2f(eigen_vecs[1].x * eigen_val[1], eigen_vecs[1].y * eigen_val[1]);
	//to map angle between -pi/2 to pi/2
	ellipse.angle = atan2(eigen_vecs[0].y, eigen_vecs[0].x); // orientation in radians
	if (ellipse.angle < 0) {
		ellipse.angle = pi - std::abs(ellipse.angle);
	}

	ellipse.longR = length(p1 - ellipse.center);
	ellipse.shortR = length(p2 - ellipse.center);

}

#if 0
/* Given a ellipse, find the minimum area ellipse that covers 95% of fiber centers (search around the given ellipse) */
void CrossSection::minAreaEllipse(const yarnIntersect2D &pts, const Ellipse &ellipse, Ellipse &minEllipse) {
	//Find the total number of points for all plys
	int sz = 0;
	for (int p = 0; p < pts.size(); ++p)
		sz += pts[p].size();
	std::vector<vec2f> data_pts (sz) ;

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

void CrossSection::extractCompressParam(const std::vector<yarnIntersect2D> &allPlaneIntersect, std::vector<Ellipse> &ellipses) {
	
	for (int i = 0; i < allPlaneIntersect.size(); ++i)
	{
		Ellipse ell;
		fitEllipse(allPlaneIntersect[i], ell);

		//Ellipse minEll;
		//minAreaEllipse(allPlaneIntersect[i], ell, minEll);
		
		ellipses.push_back(ell);
	}

	std::cout << "Compression parameters for each cross-sections are written to the file! \n";
}


void CrossSection::parameterizeEllipses(const std::vector<Ellipse> &ellipses, std::vector<Ellipse> &simple_ellipses) {

	const int ignorPlanes = 0.15 * ellipses.size();
	float ell_lng_avg = 0.f, ell_shrt_avg = 0.f, ell_angle_avg = 0.f;
	for (int i = ignorPlanes; i < ellipses.size() - ignorPlanes; ++i)
	{
		ell_lng_avg += ellipses[i].longR;
		ell_shrt_avg += ellipses[i].shortR;
		ell_angle_avg += std::abs(ellipses[i].angle);
	}
	ell_lng_avg /= (ellipses.size() - 2 * ignorPlanes);
	ell_shrt_avg /= (ellipses.size() - 2 * ignorPlanes);
	ell_angle_avg /= (ellipses.size() - 2 * ignorPlanes);



	//unsigned int numberOfPoints = ellipses.size();
	//Point2DVector points;
	//for (int i = ignorPlanes; i < ellipses.size() - ignorPlanes; ++i) {
	//	Eigen::Vector2d point;
	//	point(0) = i;
	//	point(1) = ellipses[i].shortR;
	//	points.push_back(point);
	//}
	//Eigen::VectorXd x(1);
	//x.fill(0.05f);
	//MyFunctorNumericalDiff functor;
	//functor.Points = points;
	//Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);
	//Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(x);


	for (int i = 0; i < ellipses.size(); ++i) {
		Ellipse ell;
		ell.longR = ell_lng_avg;
		ell.shortR = ell_shrt_avg;
		if (i > ignorPlanes && i < ellipses.size() - ignorPlanes)
			ell.angle = ell_angle_avg;
		else
			ell.angle = ellipses[i].angle; //since we don't care about the two ends
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

	/******************************************************************/
	world_point.y = dot(vec3f(e1.x, e2.x, n.x), local) + plane.point.x;
	world_point.x = dot(vec3f(e1.y, e2.y, n.y), local) + plane.point.y;
	world_point.z = dot(vec3f(e1.z, e2.z, n.z), local) + plane.point.z;
	//note that order of the axis matches with the project2plane()
}


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

void CrossSection::deCompressYarn(const std::vector<yarnIntersect2D> &planeIts, std::vector<Ellipse> &ellipses, std::vector<yarnIntersect2D> &deCompressPlaneIts) {
	deCompressPlaneIts.resize(planeIts.size() );
	for (int i = 0; i < planeIts.size(); ++i) // number of planes
	{
		deCompressPlaneIts[i].resize(planeIts[i].size() );
		for (int p = 0; p < planeIts[i].size(); ++p) // number of plys
		{
			deCompressPlaneIts[i][p].resize(planeIts[i][p].size() );
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
				_p_long *= 0.0286676 / ellipse_long;
				_p_short *= 0.0286676 / ellipse_short;
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

	for (int i = 0; i < e1e2_Its.size(); ++i) { // all planes
		vec3f n = m_planesList[i].n;
		vec3f e1 = m_planesList[i].e1;
		vec3f e2 = m_planesList[i].e2;
		int ply_num = e1e2_Its[i].size();
		for (int p = 0; p < ply_num; ++p) { // all plys
			vec2f plyCntr(0.f);
			for (int v = 0; v < e1e2_Its[i][p].size(); ++v) { // ply intersections
				vec2f e1e2_p= e1e2_Its[i][p][v];

				/* transform plyCenters in e1-e2 coord to xy plane */
				/****************************************************/
				// project e2 to ex, and e1 to ey with no translation
				xy_Its[i][p][v].y = dot(vec2f(e1.x, e2.x), e1e2_p);
				xy_Its[i][p][v].x = dot(vec2f(e1.y, e2.y), e1e2_p);
			}
		}
	}
}

void CrossSection::extractPlyTwist(const std::vector<yarnIntersect2D> &decompressPlaneIts, const char *plyCenterFile) {

	std::ofstream fout(plyCenterFile);
	for (int i = 0; i < decompressPlaneIts.size(); ++i ) { // all planes	
		//first find the yarn center
		//No need to yarnCenter because planes are generated using their center point
		//vec2f yarnCntr = vec2f(m_planesList[i].point.x, m_planesList[i].point.y);

		int ply_num = decompressPlaneIts[i].size();
		for (int p=0; p<ply_num; ++p ) { // all plys
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
	std::ofstream fout(ParameterizePlyCntrFile);
	assert(!fout.fail());

	for (int i = 0; i < m_planesList.size(); ++i) {

		vec2f plyCntr1(0.f), plyCntr2(0.f);
		fin >> plyCntr1.x >> plyCntr1.y;
		fin >> plyCntr2.x >> plyCntr2.y;

		//find R between yarnCenter and plyCenter[0]
		float R = length(plyCntr1);
		float theta = atan2(plyCntr1.y, plyCntr1.x);
		allR.push_back(R);
		allTheta.push_back(theta);	
	}

	// use LM to fit a curve to R and theta
	const int ignorPlanes = 0.15 *  m_planesList.size();
	unsigned int numberOfPoints = m_planesList.size();
	Point2DVector points;
	float R_avg = 0.f;
	for (int i = ignorPlanes; i < m_planesList.size() - ignorPlanes; ++i) {
		Eigen::Vector2d point;
		point(0) = i;
		point(1) = allR[i];
		points.push_back(point);
		R_avg += allR[i];
	}
	R_avg /= static_cast<float>(m_planesList.size() - 2 * ignorPlanes);
	Eigen::VectorXd x(4);
	//x.fill(0.05f);
	x(0) = 0.01;
	x(1) = pi / 25.f;
	x(2) = 0.f;
	x(3) = 0.01;

	MyFunctorNumericalDiff functor;
	functor.Points = points;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);
	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(x);

	fout << m_planesList.size() << '\n';
	for (int i = 0; i < allR.size(); ++i) {
		float R = x(0) * sin(x(1)*i + x(2)) + x(3);
		//float R = 0.01 * sin((pi/25.f) * i + 0) + 0.01;
		fout << R_avg << " " << allTheta[i] << '\n';
	}
	fout.close();
}

void CrossSection::extractFiberVectors(const std::vector<yarnIntersect2D> &decompressPlaneIts, std::vector<plyItersect2D> &fiberCntrVector) {
	fiberCntrVector.resize(decompressPlaneIts.size());
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
	for (int i = 0; i < fiberCntrVector.size() - 1; ++i) { // all planes (-1 because it checks the angle between two consecutive points)
		float theta_avg = 0;
		const int fiber_num = std::min( fiberCntrVector[i].size(), fiberCntrVector[i+1].size() ); //At the endpoints, next plane might not have as many as the current plane
		for (int v = 0; v < fiber_num; ++v) { // fiber-centers for each plane
			//compute the angle between two vec2f
			float cos = dot(fiberCntrVector[i][v], fiberCntrVector[i + 1][v]) / (length(fiberCntrVector[i + 1][v]) * length(fiberCntrVector[i][v]));
			if (1.f - cos > 1e-7) // skip the nan that are generated duo to small acos()
				theta_avg += acos(cos);
		}
		
		theta_avg /= fiber_num;
		fiber_theta.push_back(theta_avg);
		fout << fiber_theta[i+1] << '\n';
	}
	fout.close();
}
