#pragma once
#ifndef _CROSS_SECTION_H_
#define _CROSS_SECTION_H_

#include "hermiteCurve.h"
#include "Fiber.h"

class CrossSection {
public:

	CrossSection(const char* yarnfile, const char* curvefile, const int subdiv) {
		init(yarnfile, curvefile, subdiv);
	}
	void init(const char* yarnfile, const char* curvefile, const int subdiv);
	void buildPlanes(std::vector<Eigen::Vector4d> planesList);
	void findIntersection(const Eigen::Vector4d &plane, std::vector<Eigen::Vector3d> pnts); //plane vector is [a,b,c,d] for ax+by+cz+d=0
	//void findIntersections(const Eigen::Vector4d &plane, std::vector<Eigen::Vector3d> pnts);

protected:
	HermiteCurve m_curve;
	Fiber::Yarn m_yarn;


};

#endif // !_CROSS_SECTION_H_
