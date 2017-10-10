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

	CrossSection(const char* yarnfile, const char* curvefile, const int subdiv) {
		init(yarnfile, curvefile, subdiv);
	}
	void init(const char* yarnfile, const char* curvefile, const int subdiv);
	void buildPlanes();
	bool linePlaneIntersection(const vec3f &start, const vec3f &end, const Plane &plane, vec3f &its);
	bool yarnPlaneIntersection(const Plane &plane, std::vector<vec3f> &its);

protected:
	HermiteCurve m_curve;
	Fiber::Yarn m_yarn;
	std::vector<Eigen::Vector4d> planesList;

};

#endif // !_CROSS_SECTION_H_
