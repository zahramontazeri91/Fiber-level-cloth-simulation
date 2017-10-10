#include "CrossSection.h"

void CrossSection::init(const char* yarnfile, const char* curvefile, const int subdiv) {
	m_yarn.build(yarnfile);
	m_curve.init(curvefile, subdiv);
	buildPlanes();

}

void CrossSection::buildPlanes() {

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
	std::cout << "Num intersection with the plane: " << itsList.size() << std::endl;
	if (isIntrsct)
		return true;
	return false;
}