#pragma once
#ifndef _CROSS_SECTION_H_
#define _CROSS_SECTION_H_

#include "hermiteCurve.h"
#include "Fiber.h"
struct Plane {
	vec3f normal;
	vec3f point;
};

class CrossSection {
public:

	CrossSection(const char* yarnfile, const int ply_num, const char* curvefile, const int subdiv, const int num_planes) {
		init(yarnfile, ply_num, curvefile, subdiv, num_planes);
	}
	void init (const char* yarnfile, const int ply_num, const char* curvefile, const int subdiv, const int num_planes);
	void buildPlanes (const int num_planes);
	bool linePlaneIntersection (const vec3f &start, const vec3f &end, const Plane &plane, vec3f &its);
	bool yarnPlaneIntersection (const Plane &plane, std::vector<vec3f> &itsList);
	bool allPlanesIntersections (std::vector<std::vector<vec3f>> &itsLists);
	void write_PlanesIntersections(const char* filename, std::vector<std::vector<vec3f>> &itsLists);
	inline void get_plane(const int i, Plane &plane) {
		plane = m_planesList[i];
	}

protected:
	HermiteCurve m_curve;
	Fiber::Yarn m_yarn;
	std::vector<Plane> m_planesList;
};

#endif // !_CROSS_SECTION_H_
