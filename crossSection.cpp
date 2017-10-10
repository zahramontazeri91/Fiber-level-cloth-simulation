#include "CrossSection.h"

void CrossSection::init(const char* yarnfile, const char* curvefile, const int subdiv) {
	m_yarn.build(yarnfile);
	m_curve.init(curvefile, subdiv);
}

void CrossSection::buildPlanes(std::vector<Eigen::Vector4d> planesList) {

}

void CrossSection::findIntersection(const Eigen::Vector4d &plane, std::vector<Eigen::Vector3d> pnts) { //TODO: change this to 2d points?

}