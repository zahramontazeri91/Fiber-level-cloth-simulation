#include "CrossSection.h"

CrossSection::CrossSection(const Fiber::Yarn &yarn) {
	m_yarn = yarn;
	m_planesList.resize(yarn.getStepNum());
	for (int step_id = 0; step_id < m_planesList.size(); step_id++) {
		const float z = yarn.getStepSize() * (step_id - yarn.getStepNum() / 2.f); // devided by 2 Bcuz yarn lies between neg and pos z
		m_planesList[step_id].point = vec3f(0.f, 0.f, z);
		m_planesList[step_id].normal = vec3f(0.f, 0.f, 1.f); //always parallel to xy plane
		m_planesList[step_id].binormal = vec3f(1.f, 0.f, 0.f);
	}
	// no need to initialize m_curve since we already have the cross-section planes
}

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
	for (int i = 0; i < num_planes; ++i) {
		Eigen::Vector3d curve_p = m_curve.eval(i * crossSection_t);
		Eigen::Vector3d curve_t = m_curve.evalTangent(i * crossSection_t);
		Eigen::Vector3d curve_n = m_curve.evalNormal(i * crossSection_t);
		Eigen::Vector3d curve_b = curve_t.cross(curve_n); //need binormal of the plane to project 3D point on the plane
		Plane plane;
		plane.point = vec3f(curve_p[0], curve_p[1], curve_p[2]);
		plane.normal = vec3f(curve_t[0], curve_t[1], curve_t[2]);
		plane.binormal = vec3f(curve_b[0], curve_b[1], curve_b[2]);
		m_planesList.push_back(plane);
		//std::cout << "plane point: " << plane.point.x << "  " << plane.point.y << "  " << plane.point.z << std::endl;
	}
}


bool CrossSection::linePlaneIntersection(const vec3f &start, const vec3f &end, const Plane &plane, vec3f &its) {
	vec3f dir = end - start;
	if (!dot(dir, (plane.normal)))
		return false;

	const double t = dot(plane.normal, (plane.point - start)) / dot(plane.normal, dir);
	vec3f hit = start + t*dir;
	// return only if it's within the segment
	if (t <= 1.0 && t >= 0.0) {
		its = hit;
		return true;
	}
	return false;
}

bool CrossSection::yarnPlaneIntersection(const Plane &plane, yarnIntersect &itsList) {
	bool isIntrsct = false;
	const int ply_num = m_yarn.plys.size();
	itsList.resize(ply_num);

	assert( nv::length(plane.normal) - 1.f < EPS   && "normal vector is not normalized!" );

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
			fprintf_s(fout, "center: %.4lf %.4lf %.4lf \n", m_planesList[cs].point[0], m_planesList[cs].point[1], m_planesList[cs].point[2]);
			/* Then all the intersections for each ply */
			for (int p = 0; p < ply_num; ++p) { //for each ply 
				fprintf_s(fout, "ply_fiber_num: %d \n", itsLists[cs][p].size());
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
	vec3f n = plane.normal;
	vec3f e1 = plane.binormal;
	vec3f e2 = cross(n, e1);

	assert(nv::length(n) - 1.f < EPS   && "Normal vector is not normalized!");
	assert(dot(n, e1) < EPS);

	P2d.x = dot(e1, (P3d -  plane.point));
	P2d.y = dot(e2, (P3d - plane.point));
	
}

void CrossSection::PlanesIntersections2D(std::vector<yarnIntersect> &itsLists, std::vector<yarnIntersect2D> &allPlaneIntersect) {
	// TODO: no need to this text file. remove it
	if (m_planesList.size() != itsLists.size())
		std::cout << itsLists.size() << " out of " << m_planesList.size() << " many planes had intersections! \n";

	const int ply_num = m_yarn.plys.size();

	for (int cs = 0; cs < itsLists.size(); ++cs) { //for each plane
		Plane plane;
		get_plane(cs, plane);
		/* First write the yarn-center */
		vec2f center2D;
		project2Plane(plane.point, plane, center2D);

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

void CrossSection::fitCircle(const yarnIntersect2D &pts, float &radius)
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
	vec2f center = vec2f(pca_analysis.mean.at<float>(0, 0), pca_analysis.mean.at<float>(0, 1));
	//Store the eigenvalues and eigenvectors
	std::vector<vec2f> eigen_vecs(1);
	//std::vector<float> eigen_val(1);
	eigen_vecs[0] = vec2f(pca_analysis.eigenvectors.at<float>(0, 0),
			pca_analysis.eigenvectors.at<float>(0, 1));

	//find eigen values by projecting the points on eigenVectors
	float max = std::numeric_limits<float>::min();
	for (int i = 0; i < data_pts.rows; ++i) {
		float prj = eigen_vecs[0].x * data_pts.at<float>(i, 0) + eigen_vecs[0].y * data_pts.at<float>(i, 1);
		if (std::abs(prj) > max)
			max = std::abs(prj);
	}
	radius = max;
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
	eigen_val[0] = std::max(max0,max1); //TODO
	eigen_val[1] = std::min(max1,max0);

	vec2f p1 = ellipse.center + vec2f(eigen_vecs[0].x * eigen_val[0], eigen_vecs[0].y * eigen_val[0]);
	vec2f p2 = ellipse.center - vec2f(eigen_vecs[1].x * eigen_val[1], eigen_vecs[1].y * eigen_val[1]);
	ellipse.angle = atan2(eigen_vecs[0].y, eigen_vecs[0].x); // orientation in radians
	ellipse.longR = length(p1 - ellipse.center);
	ellipse.shortR = length(p2 - ellipse.center);
}

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
void CrossSection::extractCompressParam(const std::vector<yarnIntersect2D> &allPlaneIntersect, std::vector<Ellipse> &ellipses, const char* compressFile) {

	for (int i = 0; i < allPlaneIntersect.size(); ++i)
	{
		Ellipse ell;
		fitEllipse(allPlaneIntersect[i], ell);
		
		Ellipse minEll;
		minAreaEllipse(allPlaneIntersect[i], ell, minEll);
		
		ellipses.push_back(ell);
	}

	//write to compress.txt (a, b , alpha)

	FILE * fout;
	if (fopen_s(&fout, compressFile, "wt") == 0) {
		fprintf_s(fout, "%d \n", allPlaneIntersect.size());
		for (int i = 0; i < allPlaneIntersect.size(); ++i)
		{
			fprintf_s(fout, "%.4f %.4f %.4f \n", ellipses[i].longR, ellipses[i].shortR, ellipses[i].angle);
		}
		fclose(fout);
	}

	std::cout << "Compression parameters for each cross-sections are written to the file! \n";
}

void CrossSection::planeIts2world(Plane &plane, vec2f &plane_point, vec3f &world_point) {
	//Finding 3D coord of a point of the plane having its 2D: inverse of project2plane()
	vec3f n = plane.normal;
	vec3f e1 = plane.binormal;
	vec3f e2 = cross(n, e1);

	vec3f local(plane_point.x, plane_point.y, 0.f); //as the point exists on the plane
	//world = [e1 e2 n] * local
	//local = [e1; e2; n] * world
	world_point.x = dot(vec3f(e1.x, e2.x, n.x), local) + plane.point.x;
	world_point.y = dot(vec3f(e1.y, e2.y, n.y), local) + plane.point.y;
	world_point.z = dot(vec3f(e1.z, e2.z, n.z), local) + plane.point.z;
	//note that order of the axis matches with the project2plane()

}

void CrossSection::extractNormals(std::vector<Ellipse> &ellipses, std::vector<vec3f> &normals) {
	for (int i = 0; i < m_planesList.size(); ++i) {
		//rotate plane.e1=[1,0] by theta in the plane to obtain ellipse-short-axis
		vec2f p2D(cos(ellipses[i].angle), sin(ellipses[i].angle));
		//now project it to world coord
		vec3f end3D;
		planeIts2world(m_planesList[i], p2D, end3D);
		vec3f start3D = m_planesList[i].point;
		//note that ellipse center is same as plane.point
		vec3f normal = end3D - start3D;
		assert(nv::length(normal) - 1.f < EPS   && "Normal vector is not normalized!");
		normals.push_back(normal);
	}
	/*
	//test planeIts2world
	Plane plane;
	plane.point = vec3f(1.f,1.f,5.f);
	plane.normal = vec3f(0, 0, 1.f);
	vec2f p2(1.f, 1.f);
	vec3f p3;
	planeIts2world(plane, p2, p3);
	std::cout << p3.x << "  " << p3.y << "  " << p3.z << std::endl;*/
}

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
				
				const float theta = ellipses[i].angle;

				// rotate points by neg theta
				vec2f rot_axis_x(1.f, 0.f), rot_axis_y(0.f, 1.f);   //rotated the axis by -theta
				//vec2f rot_axis_x = rot2D(axis_x, -1.f * theta);
				//vec2f rot_axis_y = rot2D(axis_y, -1.f * theta);
				//assert(nv::dot(rot_axis_x, rot_axis_y) == 0);
				float _p_x = nv::dot(vec2f(its.x, its.y), rot_axis_x);
				float _p_y = nv::dot(vec2f(its.x, its.y), rot_axis_y);

				// scale the shape of cross section
				// TODO: scale to yarn.radius()
				_p_x *= 0.0286676 / ellipse_long;
				_p_y *= 0.0286676 / ellipse_short;

				vec2f new_p = _p_x * rot_axis_x + _p_y * rot_axis_y;

				deCompressPlaneIts[i][p][v] = vec2f(new_p.x, new_p.y);
			}
		}
	}
}

void CrossSection::extractPlyTwist(const std::vector<yarnIntersect2D> &allPlaneIntersect, const char *plyCenterFile) {
	//std::vector<std::vector<float>> helixRad;
	//std::vector<std::vector<float>> helixTheta;
	//helixRad.resize(allPlaneIntersect.size());
	//helixTheta.resize(allPlaneIntersect.size());

	std::ofstream fout(plyCenterFile);
	for (int i = 0; i < allPlaneIntersect.size(); ++i ) { // all planes	
		//first find the yarn center
		//No need to yarnCenter because planes are generated using their center point
		//vec2f yarnCntr = vec2f(m_planesList[i].point.x, m_planesList[i].point.y);

		int ply_num = allPlaneIntersect[i].size();
		for (int p=0; p<ply_num; ++p ) { // all plys
			vec2f plyCntr(0.f);
			for (int v = 0; v < allPlaneIntersect[i][p].size(); ++v) { // ply intersections
				plyCntr += allPlaneIntersect[i][p][v];
			}
			plyCntr /= static_cast<float>(allPlaneIntersect[i][p].size());
			//vec2f plyVec = plyCntr - yarnCntr;
			
			fout << plyCntr.x << " " << plyCntr.y << '\n';
			
			//float radius = std::sqrt ( std::pow(plyVec.x, 2.0) + std::pow(plyVec.y, 2.0) );
			//float theta = atan2(plyVec.y, plyVec.x) - 2 * pi * static_cast<float>(i) / allPlaneIntersect.size();
			//helixRad[i].push_back(radius);
			//helixTheta[i].push_back(theta);
		}
		fout << '\n';
	}
	fout.close();
}

