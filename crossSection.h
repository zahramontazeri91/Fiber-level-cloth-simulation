#pragma once
#ifndef _CROSS_SECTION_H_
#define _CROSS_SECTION_H_

#include <opencv2/opencv.hpp> 
#include "hermiteCurve.h"
#include "Fiber.h"

#define epsilon std::numeric_limits<float>::epsilon()

typedef std::vector<vec3f> plyItersect;				     //Plane intersection with each of the plys
typedef std::vector<plyItersect> yarnIntersect;		     //Plane intersection with whole yarn
typedef std::vector<vec2f> plyItersect2D;				     //Plane intersection with each of the plys in 2D
typedef std::vector<plyItersect2D> yarnIntersect2D;		     //Plane intersection with whole yarn in 2D


struct Plane {
	Plane() : normal(vec3f(0.f)), binormal(vec3f(0.f)), point(vec3f(0.f)) {}
	vec3f normal;    //spline tangent
	vec3f binormal;  //spline normal  //TODO: check if this should be spline binormal!
	vec3f point;
};


struct Ellipse {
	Ellipse() : center(vec2f(0.f)), longR(0.f), shortR(0.f) {}
	vec2f center;
	float longR;
	float shortR;
	float angle;
};

class CrossSection {
public:

	CrossSection(const char* yarnfile, const int ply_num, const char* curvefile, const int subdiv, const int num_planes) {
		init(yarnfile, ply_num, curvefile, subdiv, num_planes);
	}
	void init (const char* yarnfile, const int ply_num, const char* curvefile, const int subdiv, const int num_planes);
	void buildPlanes (const int num_planes);
	/* Intersection between a segment, defined between start to end, with a plane */
	bool linePlaneIntersection (const vec3f &start, const vec3f &end, const Plane &plane, vec3f &its);
	/* Intersection between the yarn (multiple plys) and a plane, Results will be stored in a vector< vector<vec3f> >  */
	bool yarnPlaneIntersection (const Plane &plane, yarnIntersect &itsList);
	/* Intersection between the yarn and all of the planes */
	bool allPlanesIntersections (std::vector<yarnIntersect> &itsLists);
	/* Write 3d points of the cross-sections with all planes to a file */
	void write_PlanesIntersections3D(const char* filename, std::vector<yarnIntersect> &itsLists);
	/* Find the 2D cross-section points of a yarn with all the planes*/
	void PlanesIntersections2D(std::vector<yarnIntersect> &itsLists, std::vector<yarnIntersect2D> &allPlaneIntersect);
	/* Project a 3D point to a 2D plane */
	void project2Plane(const vec3f& P3d, const Plane& plane, vec2f& P2d);
	/* Return the ith plane */
	inline void get_plane(const int i, Plane &plane) {
		plane = m_planesList[i];
	}
	/* For 2D points gathered as an ellipse, return eigen values and eigen vectors in ellipse format */
	void getOrientation(const yarnIntersect2D &pts, Ellipse &ellipse);
	/* Get ellipse a, b and angle for each cross-section and write it to the file */
	void extractCompressParam(const std::vector<yarnIntersect2D> &allPlaneIntersect, std::vector<Ellipse> &ellipses, const char* filename);
	/* Given a ellipse, find the minimum area ellipse that covers 95% of fiber centers (search around the given ellipse) */
	void CrossSection::minAreaEllipse(const yarnIntersect2D &pts, const Ellipse &ellipse, Ellipse &minEllipse);

protected:
	HermiteCurve m_curve;
	Fiber::Yarn m_yarn;
	std::vector<Plane> m_planesList;
};

#endif // !_CROSS_SECTION_H_
