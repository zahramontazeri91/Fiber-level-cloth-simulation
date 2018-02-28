#pragma once
#include "crossSection.h"
#include "Fiber.h"

//void interpolate(const std::vector<float> vals, const int interval, std::vector<float> newVals);
//void appendCompress_yarn(const std::vector<Fiber::Yarn::Compress> &compress_segs, const int seg_vrtx, const int yarn_vrtx, const char* compressFile);
//void appendCenter_yarn(const std::vector<Fiber::Yarn::CenterLine> &centerlines, const int seg_vrtx, const int yarn_vrtx, const char* curveFile);

void decomposeS(const Matrix_S &mat_S, Ellipse &ellipse);

void writeParameters(std::vector<Matrix_S> &all_mat_S, std::vector<float> &all_theta_R, const char* compress_R, const char* compress_S);

void extractCompress_seg(const char* configfile, const char* yarnfile1, const char* yarnfile2, const char* deformGrad, const char* compress_S,
	const char* curveFile, const char* normFile, const int ply_num, const int vrtx_num );

//void constFitting_compParam(const std::vector<Ellipse> &ellipses, const std::vector<float> &theta_R,
	//const float trimPercent, Fiber::Yarn::Compress &compress);
//void sinFitting_curve(const char* curveFile, const float trimPercent, Fiber::Yarn::CenterLine &curve);

//void fittingPlyCenter(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* plyCenterFile);
//void fittingCompress(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const char* compressFile);
//void fittingFiberTwisting(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* fiberTwistFile);
void plotIntersections(const std::vector<yarnIntersect2D> &its, const char* filename, const float trimPercent);
void deformRef(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Ellipse> &ellipses, const std::vector<float> &all_theta_R);
void deformRef(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Matrix_S> &all_mat_S, const std::vector<float> &all_theta_R);
void deformRef(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Eigen::Matrix2f> &all_A);

void L2norm(const std::vector<yarnIntersect2D> &its_deform, const std::vector<yarnIntersect2D> &its_trans, std::vector<float> &L2, const char* filename, const float trimPercent);

void nonPeriodicTheta(const std::vector<float> &theta, std::vector<float> &theta_new);


/****** Preparing data for training *****/
void assign_dg(const char* physical_world, std::vector<Eigen::Matrix3f> &all_world_dg);
void assign_S(const char* compress_S, std::vector<Eigen::Matrix2f> &all_S);
void transfer_dg_2local(std::vector<Eigen::Vector3d> &all_tang, std::vector<Eigen::Vector3d> &all_norm,
	std::vector<Eigen::Matrix3f> &all_world_dg, std::vector<Eigen::Matrix3f> &all_local_dg);
void rotate_S_2local (Eigen::Matrix2f &S, Eigen::Matrix2f &S_local, float &angle);
float get_angle(Eigen::Vector3d &norm1, Eigen::Vector3d &norm2, Eigen::Vector3d &tang);