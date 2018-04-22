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
void assign_twist(const char* twistFile, std::vector<float> &twists);
void transfer_dg_2local(const std::vector<Eigen::Vector3d> &all_tang, const std::vector<Eigen::Vector3d> &all_norm,
	const std::vector<Eigen::Matrix3f> &all_world_dg, std::vector<Eigen::Matrix3f> &all_local_dg, const int flip=0);
void rotate_S_2local (const Eigen::Matrix2f &S, Eigen::Matrix2f &S_local, const float &angle, const int flip=0);
float get_angle(Eigen::Vector3d &norm1, Eigen::Vector3d &norm2, Eigen::Vector3d &tang);

// old version for training data
void buildTraning_all(Fiber::Yarn &yarn, int skipFactor, int frame0, int frame1, int yarnNum, std::string &dataset, const int window_size, const float trimPercent, const int isTrain);
void buildTraining(const char* curvefile_ds, const char* normfile_ds, const char* physical_world, const char* compress_S, Fiber::Yarn &yarn, const float trimPercent,
	const int window_size, const char* curvefile_us, const char* angles, const char* physical_local_window, const char* physical_local_window_test, const char* compress_S_window,
	std::ofstream &fout_trainX_all, std::ofstream &fout_trainY_all, int isTrain);

// full pipeline:
void step0_curveSetup(const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset);
void step1_dg2local(const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset);
void step1_shapematching(const char* yarnfile1, const char* configfile, const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset);
void step2_buildTrainData(const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const int isTrain, const int window_size, const float trimPercent, const int upsample);
void step3_appendTraining(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset);
void step4_NN_output(const char* configfile, const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const int isCompress);
void full_pipeline(const char* yarnfile1, const char* configfile, const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset,
	const int isTrain, const int window_size, const float trimPercent, const int upsample = 2);