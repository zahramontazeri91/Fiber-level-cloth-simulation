#include "CrossSection.h"


void CrossSection::init(const char* yarnfile, const int ply_num, const char* curvefile, const int seg_subdiv, const int num_planes) {
	m_yarn.build(yarnfile, ply_num);
	m_curve.init(curvefile, seg_subdiv);
	buildPlanes(num_planes);
}

void CrossSection::buildPlanes(const int num_planes) {
	const double curveLen = m_curve.totalLength();
	const double crossSectionLen = curveLen / static_cast<double>(num_planes - 1);
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

	double t = dot(plane.normal, (plane.point - start)) / dot(plane.normal, dir);
	vec3f hit = start + t*dir;
	// return only if it's within the segment
	if (dot(normalize(start - hit), normalize(end - hit)) == -1) {
		its = hit;
		return true;
	}

	return false;
}

bool CrossSection::yarnPlaneIntersection(const Plane &plane, yarnIntersect &itsList) {
	bool isIntrsct = false;
	const int ply_num = m_yarn.plys.size();
	itsList.resize(ply_num);

	assert( nv::length(plane.normal) - 1.f < epsilon   && "normal vector is not normalized!" );

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
						float dist = nv::distance(its, plane.point); //distance between the fiber and the center
						if (dist < min_dist)
							closest_its = its;
						else 
							break;
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

	if (m_planesList.size() != itsLists.size())
		std::cout << itsLists.size() << " out of " << m_planesList.size()  << " many planes had intersections! \n";
	
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

	assert(nv::length(n) - 1.f < epsilon   && "Normal vector is not normalized!");
	assert(dot(n, e1) < epsilon );

	P2d.x = dot(e1, (P3d -  plane.point));
	P2d.y = dot(e2, (P3d - plane.point));
	
}

void CrossSection::write_PlanesIntersections2D(std::vector<yarnIntersect> &itsLists, std::vector<yarnIntersect2D> &allPlaneIntersect) {
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

void CrossSection::getOrientation(const yarnIntersect2D &pts, Ellipse &ellipse)
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
		if (prj0 > max0)
			max0 = prj0;
		float prj1 = eigen_vecs[1].x * data_pts.at<float>(i, 0) + eigen_vecs[1].y * data_pts.at<float>(i, 1);
		if (prj1 > max1)
			max1 = prj1;
	}
	eigen_val[0] = max0;
	eigen_val[1] = max1;

	vec2f p1 = ellipse.center + vec2f(eigen_vecs[0].x * eigen_val[0], eigen_vecs[0].y * eigen_val[0]);
	vec2f p2 = ellipse.center - vec2f(eigen_vecs[1].x * eigen_val[1], eigen_vecs[1].y * eigen_val[1]);
	ellipse.angle = atan2(eigen_vecs[0].y, eigen_vecs[0].x); // orientation in radians
	ellipse.longR = length(p1 - ellipse.center);
	ellipse.shortR = length(p2 - ellipse.center);
}

void CrossSection::minAreaEllipse(const yarnIntersect2D &pts, const Ellipse &ellipse, Ellipse &minEllips) {
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
	int insidePnt = 0;
	for (int i = 0; i < sz; ++i) {
		float x = data_pts[i].x;
		float y = data_pts[i].y;
		float t = (x - ellipse.center.x)*(x - ellipse.center.x) / (Rx*Rx) + (y - ellipse.center.y)*(y - ellipse.center.y) / (Ry*Ry);
		if (t <= 1)
			insidePnt++;
	}
	std::cout << static_cast<float> (insidePnt) / static_cast<float> (sz) << std::endl;
	
}
void CrossSection::extractCompressParam(const std::vector<yarnIntersect2D> &allPlaneIntersect, std::vector<Ellipse> &ellipses, const char* compressFile) {
	//TO DO: write compress.txt
	for (int i = 0; i < allPlaneIntersect.size(); ++i)
	{
		Ellipse ell;
		getOrientation(allPlaneIntersect[i], ell);		
		
		Ellipse minEll;
		minAreaEllipse(allPlaneIntersect[i], ell, minEll);

		ellipses.push_back(ell);
	}
	std::cout << "Compression parameters for each cross-sections are written to the file! \n";
}
