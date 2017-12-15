#pragma once
#include "crossSection.h"
#include "Fiber.h"

void interpolate(const std::vector<float> vals, const int interval, std::vector<float> newVals);
void appendCompress_yarn(const std::vector<Fiber::Yarn::Compress> &compress_segs, const int seg_vrtx, const int yarn_vrtx, const char* compressFile);
void appendCenter_yarn(const std::vector<Fiber::Yarn::CenterLine> &centerlines, const int seg_vrtx, const int yarn_vrtx, const char* curveFile);

void extractCompress_seg(const char* yarnfile1, const char* yarnfile2, const char* compressFile, const int ply_num, const int vrtx_num);
void regularize_seg();

void fittingPlyCenter(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* plyCenterFile);
void fittingCompress(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const char* compressFile);
void fittingFiberTwisting(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* fiberTwistFile);
void plotIntersections(const std::vector<yarnIntersect2D> &its, const char* filename);
void deformRef(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Ellipse> &ellipses, const std::vector<float> &all_theta_R);
void L2norm(const std::vector<yarnIntersect2D> &its_deform, const std::vector<yarnIntersect2D> &its_trans, std::vector<float> &L2, const char* filename);
