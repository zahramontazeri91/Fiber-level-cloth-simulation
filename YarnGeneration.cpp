#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"
#include <string>
#include <cstdlib>

#include <stdio.h>
#include <cmath>


void loadSamples(const char* curvefile, std::vector<Eigen::Vector3f> &pnts) {
	std::ifstream fin(curvefile);
	assert(fin.is_open() && "curvefile file wasn't found!\n");
	int pnts_num = 0;
	fin >> pnts_num;
	Eigen::Vector3f pnt;
	for (int i = 0; i < pnts_num; i++) {
		fin >> pnt[0] >> pnt[1] >> pnt[2];
		pnts.push_back(pnt);
	}

}

void findAABB(std::vector<Eigen::Vector3f> &pnts, float &min_x, float &min_y, float &min_z, float &max_x, float &max_y, float &max_z ) {

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

	//std::cout << min_x << " " << min_y << " " << min_z << " " << max_x << " " << max_y << " " << max_z;
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
			if (idx_x + 1 != resol[0] ) 
				vol[idx_x + 1][idx_y][idx_z] = 1.f;
		if ((pnts[i][0] - bottom_x) < radius) 
			if (idx_x - 1 >= 0) 
				vol[idx_x - 1][idx_y][idx_z] = 1.f;

		if ((top_y - pnts[i][1]) < radius) 
			if (idx_y + 1 != resol[1]) 
				vol[idx_x][idx_y+1][idx_z] = 1.f;
		if ((pnts[i][1] - bottom_y) < radius) 
			if (idx_y - 1 >= 0) 
				vol[idx_x][idx_y-1][idx_z] = 1.f;

		if ((top_z - pnts[i][2]) < radius) 
			if (idx_z + 1 != resol[2]) 
				vol[idx_x][idx_y][idx_z+1] = 1.f;
		if ((pnts[i][2] - bottom_z) < radius) 
			if (idx_z - 1 >= 0) 
				vol[idx_x][idx_y][idx_z-1] = 1.f;

	}
}

int writeVol() {
	std::string dataset = "single_yarn/yarn4/teeth/4_1.6";
	const int f = 17000;
	const int y = 0;
	//for loop over yarns y
	std::string tmp2 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";  ////////////////
	const char* curvefile_us = tmp2.c_str();
	const int resol_x = 30;
	const int resol_y = 30;
	const int resol_z = 100;
	const float radius = 0.1;

	std::vector<Eigen::Vector3f> pnts;
	loadSamples(curvefile_us, pnts);

	float min_x, min_y, min_z, max_x, max_y, max_z;
	findAABB(pnts, min_x, min_y, min_z, max_x, max_y, max_z);

	float minAABB[3], maxAABB[3];
	float *data;
	float *data_vol;
	//int resol[3];
	int N;

	minAABB[0] = min_x-radius, minAABB[1] = min_y-radius, minAABB[2] = min_z-radius;
	maxAABB[0] = max_x+radius, maxAABB[1] = max_y+radius, maxAABB[2] = max_z+radius;

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
	fopen_s(&fout, "testVOL.vol", "wb");
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

	return 0;
}

int main(int argc, const char **argv) {


	const int phase = atoi(argv[1]);
	assert(phase == 0 || phase == 1 || phase == 2);

	if (phase == 0 ) {
		writeVol();
		return 0;

		Fiber::Yarn yarn;
		yarn.parse(argv[2]);
		yarn.simulate_ply_shuang();
		yarn.write_plys("test_ply.txt");
		const int K = yarn.getPlyNum();
		yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
		yarn.build("test_fly.txt", K);
		yarn.write_yarn("genYarn.txt");
		return 0;
	}

	if ( argc < 4 ) {
		std::cout << "Number of argument: " << argc << std::endl;
		printf("USAGE: YarnGeneration [phase1/phase2] [configFile] [datasetFile] -w window-size=10 -s upsample=2 -t isTrain=1 -x trimPercent -k skipfactor=500 -v vrtx-num -c isCompress \n");
		printf("EXAMPLE: YarnGeneration 1 yarnTypes/yarn4/config_step2.txt yarnTypes/yarn4/datasets.txt -w 10 -s 2 -t 1 -x 0.1 -k 500 -v 300 -c 1 -z 0.02 \n");
		return 1;
	}

	const char* configfile = argv[2];
	const char* datasetfile = argv[3];
	int window_size = 10;
	int upsample = 2;
	int isTrain = 1;
	int skipFactor = 500;
	float trimPercent = 0.1;
	int isCompress = 1;
	int vrtx_ds = 300;
	float stepSize_ds = 0.02;
	

	for (int i = 4; i < argc-1; ++i) {
		std::string arg = argv[i];
		std::stringstream arg1(argv[i+1]);
		if (arg == "-w")
			arg1 >> window_size;
		if (arg == "-s")
			arg1 >> upsample;
		if (arg == "-t")
			arg1 >> isTrain;
		if (arg == "-k")
			arg1 >> skipFactor;
		if (arg == "-x")
			arg1 >> trimPercent;
		if (arg == "-v")
			arg1 >> vrtx_ds;
		if (arg == "-c")
			arg1 >> isCompress;
		if (arg == "-z")
			arg1 >> stepSize_ds;
	}
	const float stepSize = stepSize_ds / float(upsample);
	const int vrtx = vrtx_ds * upsample;

	//if the data is for training, both ends must have some trimPercent
	if (isTrain)
		assert(trimPercent > 0 && "For training, both ends must be trimmed, choose -t larger than 0");

	std::ifstream fin0(configfile);
	std::ifstream fin00(datasetfile);
	if (phase) {
		assert(fin0.is_open() && "config file wasn't found!\n");
		assert(fin00.is_open() && "datasetfile file wasn't found!\n");
	}

	switch (phase) {
		case 1: {		

			const char* yarnfile1 = "genYarn_ref.txt";
			/* This yarn is the reference yarn for shapemaching (no flyaway) */
			Fiber::Yarn yarn_ref;
			yarn_ref.parse(configfile);
			yarn_ref.setStepNum(vrtx);
			yarn_ref.yarn_simulate();
			yarn_ref.write_yarn(yarnfile1);

			std::string line;
			while (std::getline(fin00, line))
			{
				if (line == "eof") break;
				std::vector<std::string> splits = split(line, ' ');
				assert(splits.size() == 5 && "USAGE: dataset-location first-frame last-frame first-yarn last-yarn");
				std::string dataset = splits[0];
				int frame0 = atoi(splits[1].c_str()) / skipFactor;
				int frame1 = atoi(splits[2].c_str()) / skipFactor + 1;
				int yarn0 = atoi(splits[3].c_str());
				int yarn1 = atoi(splits[4].c_str());
				full_pipeline(yarnfile1, configfile, vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent, upsample, stepSize);
			}

			break;
		}
		case 2: {
			std::string line;
			while (std::getline(fin00, line))
			{
				std::cout << " *********** \n \n";
				if (line == "eof") break;
				std::vector<std::string> splits = split(line, ' ');
				std::string dataset = splits[0];
				int frame0 = atoi(splits[1].c_str()) / skipFactor;
				int frame1 = atoi(splits[2].c_str()) / skipFactor + 1;
				int yarn0 = atoi(splits[3].c_str());
				int yarn1 = atoi(splits[4].c_str());
				step4_NN_output(configfile, vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, isCompress, stepSize);
			}
			break;
		}
	}
	//	std::system("pause"); //add breakpoint instead 
	return 0;
}