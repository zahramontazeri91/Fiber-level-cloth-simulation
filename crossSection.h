#pragma once
#ifndef _CROSS_SECTION_H_
#define _CROSS_SECTION_H_

#include "hermiteCurve.h"
#include "Fiber.h"

#define epsilon std::numeric_limits<float>::epsilon()

typedef std::vector<vec3f> plyItersect;				//Plane intersection with each of the plane
typedef std::vector<plyItersect> yarnIntersect;		//Plane intersection with whole yarn


struct Plane {
	Plane() : normal(vec3f(0.f)), binormal(vec3f(0.f)), point(vec3f(0.f)) {}
	vec3f normal;    //spline tangent
	vec3f binormal;  //spline normal  //TODO: check if this should be spline binormal!
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
	bool yarnPlaneIntersection (const Plane &plane, yarnIntersect &itsList);
	bool allPlanesIntersections (std::vector<yarnIntersect> &itsLists);
	void write_PlanesIntersections3D(const char* filename, std::vector<yarnIntersect> &itsLists);
	void write_PlanesIntersections2D(const char* filename, std::vector<yarnIntersect> &itsLists);
	void project2Plane(const vec3f& P3d, const Plane& plane, vec2f& P2d);
	inline void get_plane(const int i, Plane &plane) {
		plane = m_planesList[i];
	}

protected:
	HermiteCurve m_curve;
	Fiber::Yarn m_yarn;
	std::vector<Plane> m_planesList;
};

#endif // !_CROSS_SECTION_H_
