#pragma once
#include "crossSection.h"

void extractCompressParams(const char* configfile, const char* yarnfile1, const char* yarnfile2, const char* compressFile);
void fittingPlyCenter(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* plyCenterFile);
void fittingCompress(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const char* compressFile);
void fittingFiberTwisting(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* fiberTwistFile);
void plotIntersections(const std::vector<yarnIntersect2D> &its, const char* filename);
void deformRef(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Ellipse> &ellipses, const std::vector<float> &all_theta_R);
void L2norm(const std::vector<yarnIntersect2D> &its_deform, const std::vector<yarnIntersect2D> &its_trans, std::vector<float> &L2, const char* filename);
