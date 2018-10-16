#pragma once
#ifndef _CROSS_SECTION_H_
#define _CROSS_SECTION_H_

#include <opencv2/opencv.hpp> 
#include "hermiteCurve.h"
#include "Fiber.h"

#define EPS std::numeric_limits<float>::epsilon()

struct Plane {
	Plane() : n(vec3f(0.f)), e1(vec3f(0.f)), e2(vec3f(0.f)) {}
	vec3f n;			 //spline tangent
	vec3f e1;			 //normalized vector to ply-center for the first ply
	vec3f e2;			 //perpendicular to e1 and normal
	vec3f point;         //plane center
};

struct Matrix_S {
	Matrix_S() : S00(0.f), S01(0.f), S11(0.f) {}
	float S00;
	float S01;
	float S11;
};

class CrossSection {
public:

	/* constructor for simulated yarn */
	CrossSection(const char* yarnfile, const char* curvefile, const char* normfile, int ply_num, int plane_num, 
		int subdiv_curve, std::vector<yarnIntersect2D> &allPlaneIntersect, const bool hasNorm = false) {
		if (hasNorm) //constructor if having normfile as input
			init_norm(yarnfile, ply_num, curvefile, normfile, subdiv_curve, plane_num, allPlaneIntersect);
		else
			init(yarnfile, ply_num, curvefile, normfile, subdiv_curve, plane_num, allPlaneIntersect);
	}
	CrossSection(const char* yarnfile, const char* configfile, std::vector<yarnIntersect2D> &allPlaneIntersect) {
		m_yarn.parse(configfile);
		m_yarn.build(yarnfile, m_yarn.getPlyNum());
		yarn2crossSections(allPlaneIntersect);
	}

	void init(const char* yarnfile, const int ply_num, const char* curvefile, const char* normfile, const int subdiv, 
		const int num_planes, std::vector<yarnIntersect2D> &allPlaneIntersect);
	void init_norm(const char* yarnfile, const int ply_num, const char* curvefile, const char* normfile, const int subdiv,
		const int num_planes, std::vector<yarnIntersect2D> &allPlaneIntersect);
	void buildPlanes (const int num_planes, std::vector<yarnIntersect> &itsLists);

	/* Intersection between a segment, defined between start to end, with a plane */
	bool linePlaneIntersection (const vec3f &start, const vec3f &end, const Plane &plane, vec3f &its);
	/* Intersection between the yarn (multiple plys) and a plane, Results will be stored in a vector< vector<vec3f> >  */
	bool yarnPlaneIntersection (const Plane &plane, yarnIntersect &itsList);
	/* find ply-centers by finding the averaging the ply intersections as an array of plyCenters for each cross-sectional plane */
	void allPlyCenters(std::vector<std::vector<vec3f>> &plyCenters, std::vector<yarnIntersect> &itsLists);
	/* Intersection between the yarn and all of the planes */
	bool allPlanesIntersections (std::vector<yarnIntersect> &itsLists);
	/* Find the 2D cross-section points of a yarn with all the planes*/
	void PlanesIntersections2D(std::vector<yarnIntersect> &itsLists, std::vector<yarnIntersect2D> &allPlaneIntersect);
	/* Project a 3D point to a 2D plane */
	void project2Plane (const vec3f& P3d, const Plane& plane, vec2f& P2d);

	/* Return the ith plane */
	inline void get_plane(const int i, Plane &plane) {
		plane = m_planesList[i];
	}

	/* shape matching pass transformation matrix A */
	void shapeMatch_A(const Eigen::MatrixXf &pnt_trans, const Eigen::MatrixXf &pnt_ref, Eigen::Matrix2f &A, std::ofstream &phase_fout);
	void yarnShapeMatch_A(const yarnIntersect2D &pnts_trans, const yarnIntersect2D &pnts_ref, Eigen::Matrix2f &A, std::ofstream &phase_fout);
	void yarnShapeMatches_A(const std::vector<yarnIntersect2D> &pnts_trans, const std::vector<yarnIntersect2D> &pnts_ref,
		std::vector<Eigen::Matrix2f> &all_A, std::ofstream &phase_fout);

	/* return m_curve */
	HermiteCurve get_curve() {
		return m_curve;
	}

	/* Given a yarn dataStructure, transform it to a vector of cross-sections (for debug use)*/
	void yarn2crossSections(std::vector<yarnIntersect2D> &itsLists);

protected:
	HermiteCurve m_curve;
	Fiber::Yarn m_yarn;
	std::vector<Plane> m_planesList;
};


#endif // !_CROSS_SECTION_H_
