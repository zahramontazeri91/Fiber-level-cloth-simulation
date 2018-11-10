#pragma once
#include "crossSection.h"
#include "Fiber.h"


// utilities
void upsampleVector(const std::vector<double> &values, const int upsample, std::vector<double> &values_us);
void upsampleMatrix3d(const std::vector<Eigen::Matrix3d> &DGs, const int upsample, std::vector<Eigen::Matrix3d> &DGs_up);
void writeVec2file( const std::vector<Eigen::Vector3d> values, const char* filename);
void plotIntersections(const std::vector<yarnIntersect2D> &its, const char* filename, const float trimPercent);
void L2norm(const std::vector<yarnIntersect2D> &its_deform, const std::vector<yarnIntersect2D> &its_trans, std::vector<float> &L2, const char* filename, const float trimPercent);


/****** Preparing data for training *****/
void transfer_dg_2local(const std::vector<Eigen::Vector3d> &all_tang, const std::vector<Eigen::Vector3d> &all_norm,
	const std::vector<Eigen::Matrix3d> &all_world_dg, std::vector<Eigen::Matrix3d> &all_local_dg, const int flip=0);
void rotate_S_2local (const Eigen::Matrix2f &S, Eigen::Matrix2f &S_local, const float &angle, const int flip=0);
float get_angle(Eigen::Vector3d &norm1, Eigen::Vector3d &norm2, Eigen::Vector3d &tang);

// ***************************
// ******** phase 1 **********
// *** phase 1 *** generate training and test data
void generateNNinput(const char* configfile, const int vrtx, const int skipFactor, const int frame0, const int frame1,
	const int yarn0, const int yarn1, const std::string &dataset, const int isTrain,
	const float scaleSim, const int ws_ds, const float trimPercent, const int upsample);

// phase1 - step0: parse the simulatd files: fibersim fibers to mitsuba forrmat, yarnsim to centerlines and fe files as world-DGs. Then upsample all 
void step0_parseSimulData(const char* fibersim_in, const char* yarnsim_in, const char* DG_in, const int vrtx,
	const int fiberNum, const int isTrain, const float scaleSim, const int upsample, std::vector<Eigen::Vector3d> &pnts,
	std::vector<Eigen::Vector3d> &norms, std::vector<Eigen::Vector3d> &tangs, std::vector<double> &twists,
	std::vector<Eigen::Matrix3d> &worldDGs, const char* centerfile, const char* normfile, const char* fibersimfile);

// phase1 - step1: trasfer simulated DGs which are in world space to local space using defined frames
void step1_DG2local( const std::vector<Eigen::Vector3d> &norms, const std::vector<Eigen::Vector3d> &tangs, 
	const std::vector<Eigen::Matrix3d> &worldDGs, std::vector<Eigen::Matrix3d> &localDGs);

// phase1 - step2: Apply shapematching and extract matrix-s for each cross-section (isTrain)
void step2_shapematching(const char* configfile, const char* fiberRefFile, const char* fibersimfile, const char* centerFile,
	const char* normFile, const char* globalRot, const int ply_num, const int vrtx, std::vector<Eigen::Matrix2f> &matrixS,
	std::vector<yarnIntersect2D> &pnts_ref, std::vector<yarnIntersect2D> &pnts_trans, std::vector<Eigen::Matrix2f> &all_R);

// phase1 - step3: Generate test data files and training files (if isTrain)
void step3_buildNNfiles(const int isTrain, const int ws_ds, const float trimPercent, const int sampleRate,
	const std::vector<Eigen::Vector3d> &all_pnts, const std::vector<Eigen::Vector3d> &all_norm, const std::vector<Eigen::Vector3d> &all_tang,
	std::vector<double> &twists, std::vector<Eigen::Matrix3d> &worldDGs, std::vector<Eigen::Matrix2f> &matrixS,
	const char* trainXfile, const char* trainYfile, const char* testXfile, const char* anglefile);

// phase1 - step4: append training data in one file (if isTrain)
void step4_appendTraining(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, const std::string &dataset);

// ***************************
// ******** phase 2 **********
// phase 2 - apply NN files and generate deformed fibers

void step5_applyNNoutput(const char* configfile, const int vrtx, int skipFactor, int frame0, int frame1,
	int yarn0, int yarn1, const std::string &dataset, const int isTrain, const int isCompress, const float stepSize);

// phase 2 - step2: write bounding box around centerlines needed for volume rendering
void step6_createVOL(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, const std::string &dataset, const int resol_x, const int resol_y, const int resol_z, const float radius);

// further upsample the fibers to avoid facet artifact
void step7_upsample(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, const std::string &dataset, const int isTrain, const int isCompress, const int sampleRate);

// Generate NN files such that it inputs neighbor frames to NN to handle temporal coherency 
#if 0
void temporalTrainData(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const int isTrain);
#endif