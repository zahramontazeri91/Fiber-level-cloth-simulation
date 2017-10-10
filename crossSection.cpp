#include "CrossSection.h"

void CrossSection::init(const char* yarnfile, const char* curvefile, const int seg_subdiv, const int num_planes) {
	m_yarn.build(yarnfile);
	m_curve.init(curvefile, seg_subdiv);
	buildPlanes(num_planes);
}

void CrossSection::buildPlanes(const int num_planes) {
	const double curveLen = m_curve.totalLength();
	const double crossSectionLen = curveLen / static_cast<double>(num_planes - 1);
	const double crossSection_t = m_curve.arcLengthInvApprox(crossSectionLen);
	for (int i = 0; i < num_planes; ++i) {
		Eigen::Vector3d curve_p = m_curve.eval(i * crossSection_t);
		Eigen::Vector3d curve_n = m_curve.evalTangent(i * crossSection_t);
		Plane plane;
		plane.point = vec3f(curve_p[0], curve_p[1], curve_p[2]);
		plane.normal = vec3f(curve_n[0], curve_n[1], curve_n[2]);
		m_planesList.push_back(plane);
	}
}

bool CrossSection::linePlaneIntersection(const vec3f &start, const vec3f &end, const Plane &plane, vec3f &its) {
	vec3f dir = end - start;
	if (!dot (dir,(plane.normal) ) )
		return false;

	double t = dot(plane.normal,(plane.point - start) ) / dot( plane.normal,dir) ;
	vec3f hit = start + t*dir;
	// return only if it's within the segment
	if (dot(normalize(start - hit), normalize(end - hit)) == -1) {
		its = hit;
		return true;
	}

	return false;
}

bool CrossSection::yarnPlaneIntersection(const Plane &plane, std::vector<vec3f> &itsList) {
	
	bool isIntrsct = false;
	const int num_fibers = m_yarn.plys[0].fibers.size();

	for (int f = 0; f < num_fibers; ++f) {
		const int num_vrtx = m_yarn.plys[0].fibers[f].vertices.size();
		for (int v = 0; v < num_vrtx - 1; ++v) {
			vec3f start = m_yarn.plys[0].fibers[f].vertices[v];
			vec3f end = (m_yarn.plys[0].fibers[f].vertices[v+1]);
			vec3f its(0.f);
			if (linePlaneIntersection(start, end, plane, its)) {
				isIntrsct = true;
				itsList.push_back(its);
			}
		}
	}
	if (isIntrsct)
		return true;
	return false;
}