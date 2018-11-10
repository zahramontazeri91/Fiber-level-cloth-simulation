#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"
#include "curveFitting.h"

#include <stdio.h>
#include <cmath>
#include <iomanip>
#include <sstream>

void generateNNinput(const char* configfile, const int vrtx, const int skipFactor, const int frame0, const int frame1,
	const int yarn0, const int yarn1, const std::string &dataset, const int isTrain,
	const float scaleSim, const int ws_ds, const float trimPercent, const int upsample) {

	// reference yarn can be any of existing yarnTypes with different vrtx num
	const char* fiberRefFile = "genYarnRef.txt";
	Fiber::Yarn yarnRef;
	yarnRef.parse(configfile);
	if (isTrain) {
		/* This yarn is the reference yarn for shapemaching (no flyaway) */
		yarnRef.setStepNum(vrtx);
		yarnRef.yarn_simulate();
		yarnRef.write_yarn(fiberRefFile);
	}
	const int fiberNum = yarnRef.getFiberNum();
	const int plyNum = yarnRef.getPlyNum();

	std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
	std::cout << "@@@@@@@@@@ " << dataset << " @@@@@@@@@@ \n";

	for (int f = frame0; f < frame1 + 1; f += skipFactor) {
		for (int y = yarn0; y < yarn1; y++) {

			std::cout << "@@@@@@@@@@@@@@@@ frame " << f << " - yarn " << y << " @@@@@@@@@@@@@@@@\n";

			std::stringstream frameSS;
			frameSS << std::setfill('0') << std::setw(7) << std::to_string(f);
			std::stringstream yarnSS;
			yarnSS << std::setfill('0') << std::setw(2) << std::to_string(y);

			// input files
			std::string tmp1 = "dataSets/" + dataset + "/fiber/frame_" + frameSS.str() + "fiber_" + yarnSS.str() + ".obj";
			const char* fibersim_in = tmp1.c_str();
			std::string tmp2 = "dataSets/" + dataset + "/yarn/frame_" + frameSS.str() + "fiber_" + yarnSS.str() + ".obj";
			const char* yarnsim_in = tmp2.c_str();
			std::string tmp3 = "dataSets/" + dataset + "/yarn/frame_" + frameSS.str() + "fiber_" + yarnSS.str() + ".fe";
			const char* DG_in = tmp3.c_str();

			// output files
			std::string tmp4 = "fibersim/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* fibersimfile = tmp4.c_str();
			std::string tmp5 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* normfile = tmp5.c_str();
			std::string tmp6 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* centerfile = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/globalRot_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* globalRot = tmp7.c_str();

			// NN output files
			std::string tmp8 = "input/" + dataset + "/trainX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* trainXfile = tmp8.c_str();
			std::string tmp9 = "input/" + dataset + "/trainY_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* trainYfile = tmp9.c_str();
			std::string tmp10 = "input/" + dataset + "/testX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* testXfile = tmp10.c_str();
			std::string tmp11 = "input/" + dataset + "/angle_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* anglefile = tmp11.c_str();

			// files needed for NN loss function
			std::string tmp12 = "input/" + dataset + "/simul2D_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* simul2Dfile = tmp12.c_str();
			std::string tmp13 = "input/" + dataset + "/ref2D_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* ref2Dfile = tmp13.c_str();

			std::vector<Eigen::Vector3d> pnts, norms, tangs;
			std::vector<Eigen::Matrix3d> worldDGs, localDGs;
			std::vector<double> twists;
			const int sampleRate = upsample; // downsample the NN testfiles back to simulated resolution

			// step 0: parse simulated data
			step0_parseSimulData(fibersim_in, yarnsim_in, DG_in, vrtx, fiberNum, isTrain, scaleSim, upsample,
				pnts, norms, tangs, twists, worldDGs, centerfile, normfile, fibersimfile);

			// step 1: transfer world DGs to local
			step1_DG2local(norms, tangs, worldDGs, localDGs);

			// step 2: extract matrix S
			std::vector<Eigen::Matrix2f> matrixS;
			if (isTrain) {
				std::vector<yarnIntersect2D> pnts_ref, pnts_trans;
				std::vector<Eigen::Matrix2f> all_R;
					step2_shapematching(configfile, fiberRefFile, fibersimfile, centerfile, normfile, globalRot, 
						plyNum, vrtx, matrixS, pnts_ref, pnts_trans, all_R);
					//write simulate and rotated referernce to be used for NN loss function
					writePnts(pnts_trans, all_R, simul2Dfile, 0, ws_ds, trimPercent, sampleRate);
					writePnts(pnts_ref, all_R, ref2Dfile, 1, ws_ds, trimPercent, sampleRate);
			}

			// step 3: generate NN test files and train files if isTrain
			step3_buildNNfiles(isTrain, ws_ds, trimPercent, sampleRate, pnts, norms, tangs,
				twists, worldDGs, matrixS, trainXfile, trainYfile, testXfile, anglefile);
		}
	}
	if (isTrain)
		step4_appendTraining(skipFactor, frame0, frame1, yarn0, yarn1, dataset);
}

void L2norm(const std::vector<yarnIntersect2D> &its_deform, const std::vector<yarnIntersect2D> &its_trans, std::vector<float> &L2, const char* filename, const float trimPercent) 
{
	if (its_deform.size() != its_trans.size())
		std::cout << its_deform.size() << " " << its_trans.size() << std::endl;
	assert(its_deform.size() == its_trans.size());
	FILE *fout;

	const int ignorPlanes = trimPercent * its_deform.size(); // crop the first and last 10% of the yarn
	if (fopen_s(&fout, filename, "wt") == 0) {
		fprintf_s(fout, "%d \n", its_deform.size());
		for (int i = ignorPlanes; i < its_deform.size() - ignorPlanes; ++i) { //number of planes
			float e = 0.f;
			for (int p = 0; p < its_deform[i].size(); ++p) { //number of plys
				for (int j = 0; j < its_deform[i][p].size(); ++j) { //number of intersections
					e += square_norm(its_deform[i][p][j] - its_trans[i][p][j]);
				}
			}
			L2.push_back(e);
			fprintf_s(fout, "%.6f \n", e);
		}
		fclose(fout);
	}
}

void plotIntersections(const std::vector<yarnIntersect2D> &its, const char* filename, const float trimPercent) {
	// plot cross-sections in 2D for debug purpos
	FILE *fout;
	// write the plycenters
	if (fopen_s(&fout, filename, "wt") == 0) {
		const int ignorPlanes = trimPercent * its.size(); // crop the first and last 10% of the yarn

		fprintf_s(fout, "plane_num: %d \n", its.size() - 2 * ignorPlanes);
		fprintf_s(fout, "ply_num: %d \n", its[0].size());
		fprintf_s(fout, "\n");

		for (int i = ignorPlanes; i < its.size() - ignorPlanes; ++i) { //number of planes
			fprintf_s(fout, "index_plane : %d \n", i - ignorPlanes);
			for (int p = 0; p < its[i].size(); ++p) { //number of plys
				fprintf_s(fout, "ply_fiber_num: %d \n", its[i][p].size());
				vec2f plyCenter(0.f);
				for (int j = 0; j < its[i][p].size(); ++j) { //number of intersections
					plyCenter += its[i][p][j];
				}
				plyCenter /= its[i][p].size();
				fprintf_s(fout, "plyCenter: %.4lf %.4lf \n", plyCenter.x, plyCenter.y);

				for (int j = 0; j < its[i][p].size(); ++j) { //number of intersections
					fprintf_s(fout, "%.4f %.4f \n", its[i][p][j].x, its[i][p][j].y);
				}
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
}

void transfer_dg_2local(const std::vector<Eigen::Vector3d> &all_tang, const std::vector<Eigen::Vector3d> &all_norm,
	const std::vector<Eigen::Matrix3d> &all_world_dg, std::vector<Eigen::Matrix3d> &all_local_dg, const int flip) {

	assert(all_tang.size() == all_world_dg.size() && all_norm.size() == all_world_dg.size());
	const int n = all_world_dg.size();
	all_local_dg.resize(n);
	Eigen::Vector3d ex, ey, ez;
	for (int i = 0; i < n; i++) {


		if (flip == 1) { //flip normal
			ez = all_tang[i];
			ey = -1.0*all_norm[i];
		}
		else if (flip == 2) { //flip tangent
			ez = -1.0*all_tang[i];
			ey = all_norm[i];
		}
		else if (flip == 3) { //flip tangent
			ez = -1.0*all_tang[i];
			ey = -1.0*all_norm[i];
		}
		else {
			ez = all_tang[i];
			ey = all_norm[i];
		}
		ex = ez.cross(ey);

		/** world to local **/
		Eigen::Matrix3d local_dg, M;
		M << ez[0], ez[1], ez[2],
			ey[0], ey[1], ey[2],
			ex[0], ex[1], ex[2];

		local_dg = M*all_world_dg[i];

		all_local_dg[i] = local_dg;
	}
}

void rotate_S_2local(const Eigen::Matrix2f &S, Eigen::Matrix2f &S_local, const float &angle, const int flip) {

	Eigen::Matrix2f R, A;

	R << cos(angle), -sin(angle),
		sin(angle), cos(angle);
	S_local = R*S* R.transpose();

	if (flip == 1) //flip normal
		A << -1, 0, 0, -1;
	else if (flip == 2)  //flip tangent
		A << 1, 0, 0, -1;
	else if (flip == 3)  //flip tangent and normal
		A << -1, 0, 0, -1;
	else
		A << 1, 0, 0, 1;

	S_local = A*S_local*A.transpose();

}

float get_angle(Eigen::Vector3d &norm1, Eigen::Vector3d &norm2, Eigen::Vector3d &tang) {

	// rotate the shape-matching matrix to align the new normal
	float angle = acos(norm1.dot(norm2));
	if (angle != angle) //dot product (must be 1) but might be larger than 1 and so acos return nan 
		angle = norm1.dot(norm2) > 0 ? 0.0 : pi;

	Eigen::Vector3d cross = (norm1.cross(norm2)).normalized();
	if (cross.norm() < eps) return 0.0;

	float dist = (cross - tang).norm();
	float dist_ = (cross - (-1.0*tang)).norm();
	angle = dist < dist_ ? angle : -1.0*angle;

	return angle;
}


void loadSamples(const char* curvefile, std::vector<Eigen::Vector3f> &pnts) {
	std::ifstream fin(curvefile);
	if (!fin.is_open())
		std::cout << curvefile << std::endl;
	assert(fin.is_open() && "curvefile file wasn't found!\n");
	int pnts_num = 0;
	fin >> pnts_num;
	pnts.resize(pnts_num);
	Eigen::Vector3f pnt;
	for (int i = 0; i < pnts_num; i++) {
		fin >> pnt[0] >> pnt[1] >> pnt[2];
		pnts[i] = pnt;
	}
}

void findAABB(std::vector<Eigen::Vector3f> &pnts, float &min_x, float &min_y, float &min_z, float &max_x, float &max_y, float &max_z) {

	const int sz = pnts.size();
	max_x = std::numeric_limits<float>::min();
	max_y = max_x;
	max_z = max_x;
	min_x = std::numeric_limits<float>::max();
	min_y = min_x;
	min_z = min_x;
	for (int i = 0; i < sz; i++) {
		//minimum
		if (pnts[i][0] < min_x)
			min_x = pnts[i][0];
		if (pnts[i][1] < min_y)
			min_y = pnts[i][1];
		if (pnts[i][2] < min_z)
			min_z = pnts[i][2];
		// maximum
		if (pnts[i][0] > max_x)
			max_x = pnts[i][0];
		if (pnts[i][1] > max_y)
			max_y = pnts[i][1];
		if (pnts[i][2] > max_z)
			max_z = pnts[i][2];
	}
}

void fillVolume(const std::vector<Eigen::Vector3f> &pnts, const float radius, const float minAABB[3], const float maxAABB[3], const int resol[3], std::vector<std::vector<std::vector<float>>> &vol) {

	//initialize vol
	vol.resize(resol[0]);
	for (int x = 0; x < resol[0]; x++) {
		vol[x].resize(resol[1]);
		for (int y = 0; y < resol[1]; y++) {
			vol[x][y].resize(resol[2]);
			for (int z = 0; z < resol[2]; z++) {
				vol[x][y].push_back(0.f);
			}
		}
	}

	const int sz = pnts.size();
	for (int i = 0; i < sz; i++) {
		const float len_x = maxAABB[0] - minAABB[0];
		const float len_y = maxAABB[1] - minAABB[1];
		const float len_z = maxAABB[2] - minAABB[2];

		int idx_x = ((pnts[i][0] - minAABB[0]) / len_x) * resol[0];
		int idx_y = ((pnts[i][1] - minAABB[1]) / len_y) * resol[1];
		int idx_z = ((pnts[i][2] - minAABB[2]) / len_z) * resol[2];

		if (idx_x == resol[0]) idx_x = idx_x - 1;
		if (idx_y == resol[1]) idx_y = idx_y - 1;
		if (idx_z == resol[2]) idx_z = idx_z - 1;

		vol[idx_x][idx_y][idx_z] = 1.f;

		// go d distance in all 6 directions and add neighbor voxels if needed
		float bottom_x = minAABB[0] + idx_x * (len_x / float(resol[0]));
		float top_x = bottom_x + (len_x / float(resol[0]));

		float bottom_y = minAABB[1] + idx_y * (len_y / float(resol[1]));
		float top_y = bottom_y + (len_y / float(resol[1]));

		float bottom_z = minAABB[2] + idx_z * (len_z / float(resol[2]));
		float top_z = bottom_z + (len_z / float(resol[2]));


		if ((top_x - pnts[i][0]) < radius)
			if (idx_x + 1 != resol[0])
				vol[idx_x + 1][idx_y][idx_z] = 1.f;
		if ((pnts[i][0] - bottom_x) < radius)
			if (idx_x - 1 >= 0)
				vol[idx_x - 1][idx_y][idx_z] = 1.f;

		if ((top_y - pnts[i][1]) < radius)
			if (idx_y + 1 != resol[1])
				vol[idx_x][idx_y + 1][idx_z] = 1.f;
		if ((pnts[i][1] - bottom_y) < radius)
			if (idx_y - 1 >= 0)
				vol[idx_x][idx_y - 1][idx_z] = 1.f;

		if ((top_z - pnts[i][2]) < radius)
			if (idx_z + 1 != resol[2])
				vol[idx_x][idx_y][idx_z + 1] = 1.f;
		if ((pnts[i][2] - bottom_z) < radius)
			if (idx_z - 1 >= 0)
				vol[idx_x][idx_y][idx_z - 1] = 1.f;

	}
}

void writeVol(const std::string &dataset, const int frame, const int yarn0, const int yarn1, const int resol_x, const int resol_y, const int resol_z, const float radius, const char* volumeFile) {

	std::vector<Eigen::Vector3f> pnts;
	//for loop over yarns y
	for (int y = yarn0; y < yarn1; y++) {
		std::string tmp = "input/" + dataset + "/centerYarn_" + std::to_string(frame) + "_" + std::to_string(y) + "_us.txt";
		const char* curvefile_us = tmp.c_str();
		loadSamples(curvefile_us, pnts);
	}

	float min_x, min_y, min_z, max_x, max_y, max_z;
	findAABB(pnts, min_x, min_y, min_z, max_x, max_y, max_z);


	float minAABB[3], maxAABB[3];
	float *data;
	float *data_vol;
	//int resol[3];
	int N;


	minAABB[0] = min_x - radius, minAABB[1] = min_y - radius, minAABB[2] = min_z - radius;
	maxAABB[0] = max_x + radius, maxAABB[1] = max_y + radius, maxAABB[2] = max_z + radius;


	// Modify tile scale here
	float scale = 1.0;

	int resol[3];
	resol[0] = resol_x;
	resol[1] = resol_y;
	resol[2] = resol_z;

	N = resol[0] * resol[1] * resol[2];
	data = new float[N];
	data_vol = new float[N];


	std::vector<std::vector<std::vector<float> > > volume;
	fillVolume(pnts, radius, minAABB, maxAABB, resol, volume);
	// flatten the volume
	int i = 0;
	for (int z = 0; z < resol[2]; z++) {
		for (int y = 0; y < resol[1]; y++) {
			for (int x = 0; x < resol[0]; x++) {
				data_vol[i] = volume[x][y][z];
				i++;
			}
		}
	}

	//Modidy data here
	for (int i = 0; i < N; i++) {
		data[i] = data_vol[i];
		//data[i] = 0.1;
	}


	//FILE *fout = fopen_s("testVOL.vol", "wb");
	FILE *fout;
	fopen_s(&fout, volumeFile, "wb");
	static const char tag[] = "VOL";
	fwrite(tag, 1, 3, fout);
	static const unsigned char ver = 0x3;
	fwrite(&ver, sizeof(uint8_t), 1, fout);
	int data_format = 1;
	fwrite(&data_format, sizeof(int), 1, fout);

	// Write resolution
	fwrite(resol, sizeof(int), 3, fout);

	int ch = 1;
	fwrite(&ch, sizeof(int), 1, fout);

	// Write AABB
	fwrite(minAABB, sizeof(float), 3, fout);
	fwrite(maxAABB, sizeof(float), 3, fout);

	// write voxel extent
	for (int i = 0; i < N; i++)
		fwrite(&(data[i]), sizeof(float), 1, fout);
	delete[] data;

	fclose(fout);

}

void writeFiberSim(const char* simulatedFibers, const int vrtx, const int fiberNum, const float scaledSim, const char* mitsubaFibers) {
	std::ofstream mitsubaFout(mitsubaFibers);
	std::ifstream simulFin(simulatedFibers);
	if (!simulFin.is_open()) std::cout << simulatedFibers << std::endl;
	assert(simulFin.is_open() && "simulatedFibers doesn't exist!\n");

	char v;
	float x, y, z;
	mitsubaFout << fiberNum << std::endl;
	for (int i = 0; i < fiberNum; i++) {
		mitsubaFout << vrtx << std::endl;
		for (int j = 0; j < vrtx; j++) {
			simulFin >> v >> x >> y >> z;
			mitsubaFout << x*scaledSim << " " << y*scaledSim << " " << z*scaledSim << std::endl;
		}
	}

	simulFin.close();
	mitsubaFout.close();
}

void readCenterlines(const char* simulCenterline, const int vrtx, const float scaleSim, const int upsample, std::vector<Eigen::Vector3d> &centerlines, 
	std::vector<Eigen::Vector3d> &norms, std::vector<Eigen::Vector3d> &tangs, std::vector<double> &twists) {

	std::ifstream simulFin(simulCenterline);
	if ( !simulFin.is_open() ) std::cout << simulCenterline << std::endl;
	assert(simulFin.is_open() && "simulCenterline doesn't exist!\n");

	// simulated files are downsampled
	const int vrtx_ds = int( vrtx / upsample ) ;
	std::vector<Eigen::Vector3d> pts_ds(vrtx_ds);
	std::vector<double> twist_ds(vrtx_ds);
	for (int i = 0; i < vrtx_ds; ++i) {
		char v;
		float x, y, z, t=0.0;
		std::string line;
		std::getline(simulFin, line);
		std::vector<std::string> splits = split(line, ' ');
		pts_ds[i][0] = double(std::stod(splits[1].c_str())*scaleSim);
		pts_ds[i][1] = double(std::stod(splits[2].c_str())*scaleSim);
		pts_ds[i][2] = double(std::stod(splits[3].c_str())*scaleSim);
		if (splits.size() == 5 )
			twist_ds[i] = double(std::stod(splits[4].c_str()));

		//std::cout << std::setprecision(8) << pts_ds[i][0] << " " << splits[1].c_str() <<  std::endl;
	}
	int seg_subdiv = 10;
	HermiteCurve curve;
	curve.init(pts_ds, seg_subdiv);

	// upsample the twist
	upsampleVector(twist_ds, upsample, twists);

	// assign local frames for each point
	std::vector<Eigen::Vector3d> all_pts(vrtx), all_tang(vrtx), all_norm(vrtx);
	
	curve.assign_twist(twists, all_pts, all_tang, all_norm, upsample);

	centerlines = all_pts;
	norms = all_norm;
	tangs = all_tang;
}

void readDGs(const char* DGfile, const int vrtx, const int upsample, std::vector<Eigen::Matrix3d> &worldDGs_us) {

	std::ifstream fin(DGfile);
	if (!fin.is_open()) std::cout << DGfile << " not fount! \n";
	assert(fin.is_open());

	// DGs are defined for the edges, get the average to define them for vertices
	const int vrtx_ds = int(vrtx / upsample);
	Eigen::Matrix3d DG, DG_pre;
	std::vector<Eigen::Matrix3d> worldDGs(vrtx_ds);

	for (int i = 0; i < vrtx_ds ; ++i) { 
		fin >> DG(0,0) >> DG(0, 1) >> DG(0, 2) >> DG(1, 0) >> DG(1, 1) >> DG(1, 2) >> DG(2, 0) >> DG(2, 1) >> DG(2, 2);
		worldDGs[i] = i ? 0.5*(DG + DG_pre) : DG; // because DGs are defined for edges
		DG_pre = DG;
	}

	worldDGs_us.resize(vrtx);
	upsampleMatrix3d(worldDGs, upsample, worldDGs_us);
}

void step0_parseSimulData(const char* fibersim_in, const char* yarnsim_in, const char* DG_in, const int vrtx,
	const int fiberNum, const int isTrain, const float scaleSim, const int upsample, std::vector<Eigen::Vector3d> &pnts,
	std::vector<Eigen::Vector3d> &norms, std::vector<Eigen::Vector3d> &tangs, std::vector<double> &twists, 
	std::vector<Eigen::Matrix3d> &worldDGs, const char* centerfile, const char* normfile, const char* fibersimfile) {

	std::cout << "*** step0: parse the simulated files ***\n";

	// write the fibersim file in mitsuba format, scale the fibers back
	if (isTrain)
		writeFiberSim(fibersim_in, vrtx, fiberNum, scaleSim, fibersimfile);

	// read centerline and twist, and upsample them back
	readCenterlines(yarnsim_in, vrtx, scaleSim, upsample, pnts, norms, tangs, twists);
	writeVec2file(pnts, centerfile);
	writeVec2file(norms, normfile);

	// read deformation-gradient, resample them to define for vertices instead of edges then upsample them back
	readDGs(DG_in, vrtx, upsample, worldDGs);
}


void step1_DG2local( const std::vector<Eigen::Vector3d> &norms, const std::vector<Eigen::Vector3d> &tangs, 
	const std::vector<Eigen::Matrix3d> &worldDGs, std::vector<Eigen::Matrix3d> &localDGs) {

	std::cout << "*** step1: Convert external-force to local coordinate ***\n";
	const int num_of_cores = omp_get_num_procs();

	int vrtx = norms.size();
	assert(vrtx == worldDGs.size());
	localDGs.resize(vrtx);

#pragma omp parallel for num_threads(num_of_cores)
	for (int v = 0; v < vrtx; ++v) {

		Eigen::Vector3d ez = tangs[v];
		Eigen::Vector3d ey = norms[v];
		Eigen::Vector3d ex = ez.cross(ey);

		/** local to world **/
		Eigen::Matrix3d local, world;
		world = worldDGs[v];

		Eigen::Matrix3d M;
		M << ez[0], ez[1], ez[2],
			ey[0], ey[1], ey[2],
			ex[0], ex[1], ex[2];

		local = M*world;

		//write converted parameters
		localDGs[v] = local;
	}
}

void step2_shapematching(const char* configfile, const char* fiberRefFile, const char* fibersimfile, const char* centerFile,
	const char* normFile, const char* globalRot, const int ply_num, const int vrtx, std::vector<Eigen::Matrix2f> &matrixS, 
	std::vector<yarnIntersect2D> &pnts_ref, std::vector<yarnIntersect2D> &pnts_trans, std::vector<Eigen::Matrix2f> &all_R)
{
	std::cout << "*** step2: Shapematching step to extract deformation matrix ***\n";
	// Generate non-deformed yarn
	//std::vector<yarnIntersect2D> pnts_ref;
	CrossSection cs1(fiberRefFile, configfile, vrtx, pnts_ref);

	// Generate deformed yarn
	//std::vector<yarnIntersect2D> pnts_trans;
	CrossSection cs2(fibersimfile, centerFile, normFile, ply_num, vrtx, 100, pnts_trans, true);

	assert( pnts_ref.size() == pnts_trans.size() );

	// write global rotations for phase-matching purpose
	//std::vector<Eigen::Matrix2f> all_R;
	cs2.yarnShapeMatches_A(pnts_trans, pnts_ref, matrixS, all_R);

	// Write the global rotation
	std::ofstream fout(globalRot);
	for (int l = 0; l < all_R.size(); l++) 
		fout << all_R[l](0, 0) << " " << all_R[l](0, 1) << " " << all_R[l](1, 0) << " " << all_R[l](1, 1) << std::endl;
	fout.close();

}

void storeNNData(const int isTrain, Eigen::VectorXd &dataX, Eigen::VectorXd &dataY, const int ws_us, const int sampleRate,
	const std::vector<Eigen::Vector3d> &all_tang, const std::vector<Eigen::Vector3d> &all_norm, 
	const std::vector<Eigen::Matrix3d> &all_world_dg, const Eigen::Matrix2f &S, const float &angle, const int augMode) 
{
	const int ws_ds = static_cast<int> (ws_us / sampleRate); //because we only store ds that was literally simulated
	dataX.resize(int(ws_ds * 9));
	dataY.resize(4);

	if (isTrain) {
		/************* store Y values **************/
		Eigen::Matrix2f S_local;
		rotate_S_2local(S, S_local, angle, augMode);
		dataY << S_local(0, 0), S_local(0, 1), S_local(1, 0), S_local(1, 1);
	}

	/************* store X values **************/
	std::vector<Eigen::Matrix3d> all_local_dg;
	transfer_dg_2local(all_tang, all_norm, all_world_dg, all_local_dg, augMode);
	int p = 0;
	for (int d = 0; d < ws_us; (d += sampleRate)) {
		Eigen::Matrix3d local_dg = all_local_dg[d];

		dataX(p * 9 + 0) = local_dg(0, 0);
		dataX(p * 9 + 1) = local_dg(0, 1);
		dataX(p * 9 + 2) = local_dg(0, 2);
		dataX(p * 9 + 3) = local_dg(1, 0);
		dataX(p * 9 + 4) = local_dg(1, 1);
		dataX(p * 9 + 5) = local_dg(1, 2);
		dataX(p * 9 + 6) = local_dg(2, 0);
		dataX(p * 9 + 7) = local_dg(2, 1);
		dataX(p * 9 + 8) = local_dg(2, 2);
		p++;
	}
}


void step3_buildNNfiles(const int isTrain, const int ws_ds, const float trimPercent, const int sampleRate,
	const std::vector<Eigen::Vector3d> &all_pnts, const std::vector<Eigen::Vector3d> &all_norm, const std::vector<Eigen::Vector3d> &all_tang, 
	std::vector<double> &twists, std::vector<Eigen::Matrix3d> &worldDGs, std::vector<Eigen::Matrix2f> &matrixS,
	const char* trainXfile, const char* trainYfile, const char* testXfile, const char* anglefile)
{
	std::cout << "*** step 3: Build training data  *** \n";

	const int vrtx = twists.size();
	assert(vrtx == all_norm.size() && vrtx == worldDGs.size());

	/* window-level */
	const int ignorPlanes = trimPercent * vrtx; // crop the first and last #% of the yarn
	const int ws_us = ws_ds * sampleRate; // first upsample the window-size

	const int window_num = ((vrtx - ws_us + 1) - 2 * ignorPlanes);

	std::vector<int> indices;
	for (int w = ignorPlanes; w < (vrtx - ws_us + 1) - ignorPlanes; w = w + sampleRate)
		indices.push_back(w);
	std::vector <Eigen::VectorXd> store_NN_testX_total(indices.size());
	std::vector <Eigen::VectorXd> store_NN_trainX_total(indices.size()*4); //because of augmentation *4
	std::vector <Eigen::VectorXd> store_NN_trainY_total(indices.size()*4); //because of augmentation *4
	std::vector <float> angle_total(indices.size());

	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int omp_i = 0; omp_i < static_cast<int>(indices.size()); ++omp_i) {
		const int w = indices[omp_i];

		//define a curve segment 
		const int start = w;
		const int end = w + (ws_us - 1);

		HermiteCurve segment;
		const int seg_subdiv = 10;
		segment.init_seg(all_pnts, start, end, seg_subdiv);
		std::vector<Eigen::Vector3d> all_pts_seg, all_tang_seg, all_norm_seg;
		std::vector<double>::const_iterator first = twists.begin() + start;
		std::vector<double>::const_iterator last = twists.begin() + end + 1;
		std::vector<double> twists_seg(first, last);
		segment.assign_twist(twists_seg, all_pts_seg, all_tang_seg, all_norm_seg, 1);

		std::vector<Eigen::Matrix3d> all_dg_seg(ws_us);
		int i = 0;
		for (int d = start; d <= end; d++) {
			all_dg_seg[i] = worldDGs[d];
			i++;
		}

		const int v_yarn = ceil((start + end) / 2.0);
		const int v_seg = ceil((end - start) / 2.0);
		assert(v_yarn == v_seg + start && "index for full yarn must be equal to index segment added with starting vertex");

		/************* write angle **************/
		Eigen::Vector3d norm1 = all_norm[v_yarn];
		Eigen::Vector3d norm2 = all_norm_seg[v_seg];
		Eigen::Vector3d tang = all_tang[v_yarn];
		//assert((all_tang[v_yarn] - all_tang_seg[v_seg]).norm() < eps && "tang-seg should be equal to tang-all!\n");
		float angle = get_angle(norm1, norm2, tang);
		angle_total[omp_i] = angle;

		/******** write test data *******/
		Eigen::VectorXd store_testX, store_testY; //testY is empty
		storeNNData(0, store_testX, store_testY, ws_us, sampleRate, all_tang_seg, all_norm_seg, all_dg_seg, matrixS[v_yarn], angle, 0);
		store_NN_testX_total[omp_i] = store_testX;

		// data-augmentation for training data
		if (isTrain) {
			
			Eigen::VectorXd store_trainX, store_trainY;
			storeNNData(isTrain, store_trainX, store_trainY, ws_us, sampleRate, all_tang_seg, all_norm_seg, all_dg_seg, matrixS[v_yarn], angle, 0);
			store_NN_trainX_total[omp_i * 4 + 0] = store_trainX;
			store_NN_trainY_total[omp_i * 4 + 0] = store_trainY;

			storeNNData(isTrain, store_trainX, store_trainY, ws_us, sampleRate, all_tang_seg, all_norm_seg, all_dg_seg, matrixS[v_yarn], angle, 1);
			store_NN_trainX_total[omp_i * 4 + 1] = store_trainX;
			store_NN_trainY_total[omp_i * 4 + 1] = store_trainY;

			storeNNData(isTrain, store_trainX, store_trainY, ws_us, sampleRate, all_tang_seg, all_norm_seg, all_dg_seg, matrixS[v_yarn], angle, 2);
			store_NN_trainX_total[omp_i * 4 + 2] = store_trainX;
			store_NN_trainY_total[omp_i * 4 + 2] = store_trainY;

			storeNNData(isTrain, store_trainX, store_trainY, ws_us, sampleRate, all_tang_seg, all_norm_seg, all_dg_seg, matrixS[v_yarn], angle, 3);
			store_NN_trainX_total[omp_i * 4 + 3] = store_trainX;
			store_NN_trainY_total[omp_i * 4 + 3] = store_trainY;
		}
	}

	// Now write the NN files outside the openmp loop
	// this was writing all window-size although we want the downsampled Dg's which are the actual simulated
	std::ofstream fout_testX(testXfile);
	std::ofstream fout_angle(anglefile);
	for (int l = 0; l < angle_total.size(); l++) {
		fout_angle << angle_total[l] << "\n";
	}
	for (int l = 0; l < store_NN_testX_total.size(); l++) {
		for (int p = 0; p < int(ws_ds * 9); p++) fout_testX << store_NN_testX_total[l](p) << " ";
		fout_testX << "\n";
	}
	fout_testX.close();
	fout_angle.close();

	if (isTrain) {
		std::ofstream fout_trainX(trainXfile);
		std::ofstream fout_trainY(trainYfile);
		for (int l = 0; l < store_NN_trainX_total.size(); l++) {
			for (int p = 0; p < int(ws_ds * 9); p++) fout_trainX << store_NN_trainX_total[l](p) << " ";
			fout_trainX << "\n";
		}
		for (int l = 0; l < store_NN_trainY_total.size(); l++) {
			for (int p = 0; p < 4; p++) fout_trainY << store_NN_trainY_total[l](p) << " ";
			fout_trainY << "\n";
		}
		fout_trainX.close();
		fout_trainY.close();
	}
}

void step4_appendTraining(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, const std::string &dataset) {

	std::cout << "*** step 4: Append training-data for all frames *** \n";

	std::string tmp0 = "input/" + dataset + "/trainX_all.txt";
	const char* all_trainX = tmp0.c_str();
	std::string tmp1 = "input/" + dataset + "/trainY_all.txt";
	const char* all_trainY = tmp1.c_str();
	std::ofstream fout_trainX_all(all_trainX);
	std::ofstream fout_trainY_all(all_trainY);
	std::string content_trainX = "";
	std::string content_trainY = "";

	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;

		int seg_subdiv = 10;
		for (int y = yarn0; y < yarn1; ++y) {
			std::string tmp4 = "input/" + dataset + "/trainX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* trainX = tmp4.c_str();
			std::string tmp6 = "input/" + dataset + "/trainY_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* trainY = tmp6.c_str();


			std::ifstream finX(trainX);
			assert(finX.is_open() && "trainX file wasn't found!\n");
			std::ifstream finY(trainY);
			assert(finY.is_open() && "trainY file wasn't found!\n");

			int i;
			for (i = 0; finX.eof() != true; i++) // get content of infile
				content_trainX += finX.get();
			i--;
			content_trainX.erase(content_trainX.end() - 1);     // erase last character
			finX.close();


			// trainY
			int j;
			for (j = 0; finY.eof() != true; j++) // get content of infile
				content_trainY += finY.get();
			j--;
			content_trainY.erase(content_trainY.end() - 1);     // erase last character
			finY.close();


		}
	}
	fout_trainX_all << content_trainX;                 // output
	fout_trainY_all << content_trainY;                 // output
	fout_trainX_all.close();
	fout_trainY_all.close();

}

void step5_applyNNoutput(const char* configfile, const int vrtx, int skipFactor, int frame0, int frame1,
	int yarn0, int yarn1, const std::string &dataset, const int isTrain, const int isCompress, const float stepSize) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** Testing-NN phase ***\n";
	std::cout << " @@@@@@@@@@ " << dataset << " @@@@@@@@@@ \n";
	const int num_of_cores = omp_get_num_procs();

	//// Procedural step
	/* This yarn is what will be compressed (has flyaways) */
	Fiber::Yarn yarn;
	yarn.parse(configfile);
	yarn.setStepNum(vrtx);
	yarn.setStepSize(stepSize);
	yarn.simulate_ply();
	yarn.write_plys("test_ply.txt");
	const int K = yarn.getPlyNum();
	yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
	yarn.build("test_fly.txt", K);

	const int total_frame = (frame1 - frame0) / skipFactor + 1; 
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;

		std::cout << "Generate fibers for frame <<< " << f << " >>> started... \n";

		for (int y = yarn0; y < yarn1; ++y) {

			std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* curvefile = tmp1.c_str();
			std::string tmp2 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* normfile = tmp2.c_str();
			std::string tmp6 = "input/" + dataset + "/testY_NN_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S = tmp6.c_str();

			std::string tmp11 = "input/" + dataset + "/globalRot_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* global_rot = tmp11.c_str();

			std::ifstream fin2(compress_S);
			if (isTrain)
				assert(fin2.is_open() && "testY file wasn't found!\n");
			std::ifstream fin3(curvefile);
			assert(fin3.is_open() && "curvefile file wasn't found!\n");
			std::ifstream fin4(normfile);
			assert(fin4.is_open() && "normfile file wasn't found!\n");

			///*******  write the yarn ******/
			std::string tmp3;
			if (isCompress) {
#ifndef IMPROVED_FLYAWAYS
				tmp3 = "output/" + dataset + "/genYarn_NN_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#else
				tmp3 = "output/" + dataset + "/genYarn_NN_fly_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#endif
			}
			else {
#ifndef IMPROVED_FLYAWAYS
				tmp3 = "output/" + dataset + "/genYarn_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#else
				tmp3 = "output/" + dataset + "/genYarn_fly_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#endif
			}

			const char* outfile = tmp3.c_str();
			////// Procedural step
			Fiber::Yarn yarn_compressed; //renew the yarn
			yarn_compressed = yarn;
			if (isCompress) {
				if (isTrain)
					yarn_compressed.compress_yarn_A(compress_S, global_rot);
				else 
					yarn_compressed.compress_yarn_A(compress_S);
			}
			yarn_compressed.curve_yarn(curvefile, normfile);
			yarn_compressed.write_yarn(outfile);

			//yarn_compressed.write_yarn_obj(outfile);

			/*******  Validate NN by L2-norm ******/
			//std::string tmp4 = "output/" + dataset + "/genYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			//const char* yarnfile_proc = tmp4.c_str(); //proc yarn
			//std::ifstream fin6(yarnfile_proc);
			//assert(fin6.is_open() && "yarn_proc file wasn't found!\n");
			//Fiber::Yarn yarn_proc;
			//yarn_proc.parse(configfile);
			//yarn_proc.build(yarnfile_proc, yarn_proc.getPlyNum());

			//std::string tmp5 = "fibersim/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			//const char* yarnfile_simul = tmp5.c_str();
			//std::ifstream fin5(yarnfile_simul);
			//assert(fin5.is_open() && "yarn_simul file wasn't found!\n");
			//Fiber::Yarn yarn_simul;
			//yarn_simul.parse(configfile);
			//yarn_simul.build(yarnfile_simul, yarn_simul.getPlyNum());

			//const int trimPercent = 0.15; // should match with building NN data
			//float L2;
			//yarn.L2norm_3D(yarn, yarn_proc, trimPercent, L2);
			//std::cout << "L2 is: " << L2 << std::endl;

		}
	}
}

void step6_createVOL(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, const std::string &dataset, 
	const int resol_x, const int resol_y, const int resol_z, const float radius) {

	std::cout << "*** step 6: Create volume phase ***\n";
	const int num_of_cores = omp_get_num_procs();

	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;

		std::string tmp = "output/" + dataset + "/volume_" + std::to_string(f) + ".vol";
		const char* volfile_us = tmp.c_str();
		std::cout << frame0 << " " << f << " " << volfile_us << " generation is started... \n";
		int curr_frame = f;
		writeVol(dataset, curr_frame, yarn0, yarn1, resol_x, resol_y, resol_z, radius, volfile_us);
	}
}

void upsampleYarn(const char* infile, const char* outfile, const int sampleRate) {
	std::ifstream fin(infile);
	std::ofstream fout(outfile);
	int fiberNum, vrtxNum;
	fin >> fiberNum;
	fout << fiberNum << std::endl;
	for (int i = 0; i < fiberNum; i++) {
		fin >> vrtxNum;
		const int vrtxNum_us = vrtxNum * sampleRate;
		fout << vrtxNum_us << std::endl;
		std::vector<Eigen::Vector3d> pts(vrtxNum);
		HermiteCurve curve;
		const int subdiv = 10;
		for (int i = 0; i < vrtxNum; ++i) fin >> pts[i][0] >> pts[i][1] >> pts[i][2];
		assert(pts.size() > 2);
		curve.init(pts, subdiv);
		double curveLength = curve.totalLength();


		for (int v = 0; v < vrtxNum_us; v++) {
			double len = double(v) / double(vrtxNum_us) * curveLength;
			double t = curve.arcLengthInvApprox(len);
			Eigen::Vector3d pos = curve.eval(t);
			fout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
		}
	}
	std::cout << outfile << " is written! \n";
	fout.close();
}

void step7_upsample(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, const std::string &dataset, 
	const int isTrain, const int isCompress, const int sampleRate) {

	/* given a fiber.txt file, return upsampled fiber.txt */
	std::cout << "*** step 7: Upsample fibers to avoid facet artifact *** \n";
	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;
		for (int y = yarn0; y < yarn1 ; y++) {
			std::string tmp1, tmp2;
			if (isCompress) {
				tmp1 = "output/" + dataset + "/genYarn_NN_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				tmp2 = "output/" + dataset + "/genYarn_NN_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
				const char* infile = tmp1.c_str();
				const char* outfile = tmp2.c_str();
				upsampleYarn(infile, outfile, sampleRate);
			}

			// simulated and bottomline fibers don't need upsampling
			//if (isTrain && isCompress) { //upsample the simulated yarns once
			//	tmp1 = "fibersim/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			//	tmp2 = "fibersim/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			//	infile = tmp1.c_str();
			//	outfile = tmp2.c_str();
			//	upsampleYarn(infile, outfile, sampleRate);
			//}
		}
	}
}


void upsampleVector(const std::vector<double> &values, const int upsample, std::vector<double> &values_us) {

	int N = values.size();
	values_us.resize(N*upsample);

	// Let's interpolate for upsampled vector
	int c = 0;
	for (int j = 0; j < N; j++) {
		for (int s = 0; s < upsample; s++) {
			if (j == N - 1) {
				values_us[c] = values[j];
			}
			else {
				const float w = float(s) / float(upsample);
				double value_us = (1.0 - w) * values[j] + w * values[j + 1];
				values_us[c] = value_us;
			}
			c++;
		}
	}
}

void upsampleMatrix3d(const std::vector<Eigen::Matrix3d> &DGs, const int upsample, std::vector<Eigen::Matrix3d> &DGs_up) {
	const int N = DGs.size();
	DGs_up.resize(N*upsample);
	// Let's interpolate for upsampled vector
	int c = 0;
	for (int j = 0; j < N; j++) {
		for (int s = 0; s < upsample; s++) {
			if (j == N - 1) {
				DGs_up[c] = DGs[j];
			}
			else {
				const float w = float(s) / float(upsample);
				Eigen::Matrix3d DG = (1.0 - w) * DGs[j] + w * DGs[j + 1];
				DGs_up[c] = DG;
			}
			c++;
		}
	}
}

void writeVec2file(const std::vector<Eigen::Vector3d> values, const char* filename) {
	std::ofstream fout(filename);
	fout << values.size() << std::endl;
	for (int i = 0; i < values.size(); i++) {
		fout << values[i][0] << " " << values[i][1] << " " << values[i][2] << std::endl;
	}
}


#if 0
void temporalTrainData(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const int isTrain) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** step 2: Add temporal support to training data  *** \n";
	std::cout << " @@@@@@@@@@ " << dataset << " @@@@@@@@@@ \n";

	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;
		for (int y = yarn0; y < yarn1; ++y) {

			
			// Padding for the first frame 
			int preFrame = f - skipFactor;
			int postFrame = f + skipFactor;
			if (i == 0)
				preFrame = f;
			else if (i == total_frame - 1)
				postFrame = f;

			std::string tmp0 = "input/" + dataset + "/NN/trainX_" + std::to_string(preFrame) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_pre = tmp0.c_str();
			std::string tmp1 = "input/" + dataset + "/NN/testX_" + std::to_string(preFrame) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_test_pre = tmp1.c_str();

			std::string tmp2 = "input/" + dataset + "/NN/trainX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg = tmp2.c_str();
			std::string tmp3 = "input/" + dataset + "/NN/testX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_test = tmp3.c_str();

			std::string tmp4 = "input/" + dataset + "/NN/trainX_" + std::to_string(postFrame) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_post = tmp4.c_str();
			std::string tmp5 = "input/" + dataset + "/NN/testX_" + std::to_string(postFrame) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_test_post = tmp5.c_str();

			std::string tmp6 = "input/" + dataset + "/NN/trainX_" + std::to_string(f) + "_" + std::to_string(y) + "_temporal.txt";
			const char* physical_local_seg_temporal = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/NN/testX_" + std::to_string(f) + "_" + std::to_string(y) + "_temporal.txt";
			const char* physical_local_seg_test_temporal = tmp7.c_str();


			if (isTrain) {
				std::ifstream finTrain_pre(physical_local_seg_pre);
				assert(finTrain_pre.is_open() && "trainx pre wasn't found!\n");
				std::ifstream finTrain(physical_local_seg);
				assert(finTrain.is_open() && "trainx wasn't found!\n");
				std::ifstream finTrain_post(physical_local_seg_post);
				assert(finTrain_post.is_open() && "trainx post wasn't found!\n");

				std::ofstream foutTrain(physical_local_seg_temporal);


				for (int i = 0; finTrain.eof() != true; i++) {
					std::string content_trainX = "";
					std::string content_trainX_pre = "";
					std::string content_trainX_post = "";

					std::getline(finTrain_pre, content_trainX_pre);
					std::getline(finTrain, content_trainX);
					std::getline(finTrain_post, content_trainX_post);

					foutTrain << content_trainX_pre << " " << content_trainX << " " << content_trainX_post << std::endl;

				}
				//content_trainX_temporal.erase(content_trainX_temporal.end() - 1);     // erase last character
				finTrain_pre.close();
				finTrain.close();
				finTrain_post.close();

				// write the updated training data

				foutTrain.close();
			}
			std::ifstream finTest_pre(physical_local_seg_test_pre);
			assert(finTest_pre.is_open() && "testx pre wasn't found!\n");
			std::ifstream finTest(physical_local_seg_test);
			assert(finTest.is_open() && "testx wasn't found!\n");
			std::ifstream finTest_post(physical_local_seg_test_post);
			assert(finTest_post.is_open() && "testx post wasn't found!\n");

			std::ofstream foutTest(physical_local_seg_test_temporal);


			for (int i = 0; finTest.eof() != true; i++) {
				std::string content_testX = "";
				std::string content_testX_pre = "";
				std::string content_testX_post = "";

				std::getline(finTest_pre, content_testX_pre);
				std::getline(finTest, content_testX);
				std::getline(finTest_post, content_testX_post);

				foutTest << content_testX_pre << " " << content_testX << " " << content_testX_post << std::endl ;

			}
			//content_trainX_temporal.erase(content_trainX_temporal.end() - 1);     // erase last character
			finTest_pre.close();
			finTest.close();
			finTest_post.close();

			// write the updated training data

			foutTest.close();
		}
	}
}
#endif
