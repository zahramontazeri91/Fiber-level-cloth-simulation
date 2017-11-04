#pragma once
#ifndef _CROSS_SECTION_H_
#define _CROSS_SECTION_H_

#include <opencv2/opencv.hpp> 
#include "hermiteCurve.h"
#include "Fiber.h"

#define EPS std::numeric_limits<float>::epsilon()

struct Plane {
	Plane() : n(vec3f(0.f)), e1(vec3f(0.f)), e2(vec3f(0.f)), point(vec3f(0.f)) {}
	vec3f n;		 //spline tangent
	vec3f e1;		 //long ellipse axis //TODO: check if this should be spline binormal!
	vec3f e2;		 //short ellipse axis
	vec3f point;     //plane center
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

	/* constructor for simulated yarn */
	CrossSection(const char* yarnfile, const char* curvefile, int ply_num, int plane_num, int subdiv_curve) {
		init(yarnfile, ply_num, curvefile, subdiv_curve, plane_num);
	}	

	void init (const char* yarnfile, const int ply_num, const char* curvefile, const int subdiv, const int num_planes);
	void buildPlanes (const int num_planes);
	/* Intersection between a segment, defined between start to end, with a plane */
	bool linePlaneIntersection (const vec3f &start, const vec3f &end, const Plane &plane, vec3f &its);
	/* Intersection between the yarn (multiple plys) and a plane, Results will be stored in a vector< vector<vec3f> >  */
	bool yarnPlaneIntersection (const Plane &plane, yarnIntersect &itsList);
	/* find ply-centers by finding the averaging the ply intersections as an array of plyCenters for each cross-sectional plane */
	void allPlyCenters(std::vector<std::vector<vec3f>> &plyCenters);
	/* Intersection between the yarn and all of the planes */
	bool allPlanesIntersections (std::vector<yarnIntersect> &itsLists);
	/* Write 3d points of the cross-sections with all planes to a file */
	void write_PlanesIntersections3D(const char* filename, std::vector<yarnIntersect> &itsLists);
	/* Find the 2D cross-section points of a yarn with all the planes*/
	void PlanesIntersections2D(std::vector<yarnIntersect> &itsLists, std::vector<yarnIntersect2D> &allPlaneIntersect);
	/* Project a 3D point to a 2D plane */
	void project2Plane(const vec3f& P3d, const Plane& plane, vec2f& P2d);
	/* project 2D coord back to world */
	void planeIts2world(Plane &plane, vec2f &plane_point, vec3f &world_point);
	/* Return the ith plane */
	inline void get_plane(const int i, Plane &plane) {
		plane = m_planesList[i];
	}
	/* For 2D points gathered as an ellipse, return eigen values and eigen vectors in ellipse format */
	void fitEllipse(const yarnIntersect2D &pts, Ellipse &ellipse);
	/* Get ellipse a, b and angle for each cross-section and write it to the file */
	void extractCompressParam(const std::vector<yarnIntersect2D> &allPlaneIntersect, std::vector<Ellipse> &ellipses, const char* filename);
	/* Extract normal vectors for the curve using fitted ellipse for each plane intersection */
	void extractNormals(std::vector<Ellipse> &ellipses, std::vector<vec3f> &normals, const char* pntsFile, const char* normsFile);
	/* de-compress simulated yarn */
	void deCompressYarn(const std::vector<yarnIntersect2D> &PlaneIts, std::vector<Ellipse> &ellipses, std::vector<yarnIntersect2D> &deCompressPlaneIts);
	/* get R and theta for ply-centers*/
	void extractPlyTwist(const std::vector<yarnIntersect2D> &allPlaneIntersect, const char *plyCenterFile);
	

	/* return m_curve */
	HermiteCurve get_curve() {
		return m_curve;
	}

protected:
	HermiteCurve m_curve;
	Fiber::Yarn m_yarn;
	std::vector<Plane> m_planesList;
};

#endif // !_CROSS_SECTION_H_
