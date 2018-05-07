#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"
#include <string>
#include <cstdlib>



int main(int argc, const char **argv) {


	const int phase = atoi(argv[1]);
	assert(phase == 0 || phase == 1 || phase == 2);

	if (phase == 0 ) {

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
		printf("USAGE: YarnGeneration [phase1/phase2] [configFile] [datasetFile] -w window-size=10 -s upsample=2 -t isTrain=1 -x trimPercent -k skipfactor=500 -v vrtx-num -c isCompress -z stepSize_ds -rx resolution-AABB-x -ry -rz -rad radius-AABB \n");
		printf("EXAMPLE: YarnGeneration 1 yarnTypes/yarn4/config_step2.txt yarnTypes/yarn4/datasets.txt -w 10 -s 2 -t 1 -x 0.1 -k 500 -v 300 -c 1 -z 0.02 -rx 5 -ry 5 -rz 20 -rad 0.1 \n");
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
	int resol_x = 5;
	int resol_y = 5;
	int resol_z = 20;
	float radius = 0.1;

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
		if (arg == "-rx")
			arg1 >> resol_x;
		if (arg == "-ry")
			arg1 >> resol_y;
		if (arg == "-ry")
			arg1 >> resol_y;
		if (arg == "-rad")
			arg1 >> radius;
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
				step5_createVOL( skipFactor, frame0, frame1, yarn0, yarn1, dataset, resol_x, resol_y, resol_z, radius);
			}
			break;
		}
	}
	//	std::system("pause"); //add breakpoint instead 
	return 0;
}