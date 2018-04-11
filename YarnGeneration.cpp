#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"
#include <string>
#include <cstdlib>


int main(int argc, const char **argv) {

	if (argv[0] == "-h" || argc < 4) {
		printf("USAGE: YarnGeneration [phase1/phase2] [configFile] [datasetFile] -w window-size=10 -s upsample=2 -t isTrain=1 -k skipfactor=500 -x trimPercent -v vrtx-num \n");
		printf("EXAMPLE: YarnGeneration 1 yarnTypes/yarn4/config_step2.txt yarnTypes/yarn4/datasets.txt -w 10 -s 2 -t 1 -k 500 -x 0.1 -v 300 \n");
	}
	const int phase = atoi(argv[1]);
	const char* configfile = argv[2];
	const char* datasetfile = argv[3];
	int window_size = 10;
	int upsample = 2;
	int isTrain = 1;
	int skipFactor = 500;
	float trimPercent = 0.1;
	int vrtx = 300;
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
			arg1 >> vrtx;
	}

	std::ifstream fin0(configfile);
	std::ifstream fin00(datasetfile);
	if (phase) {
		assert(fin0.is_open() && "config file wasn't found!\n");
		assert(fin00.is_open() && "datasetfile file wasn't found!\n");
	}

	switch (phase) {
		case 0: {
			const char* yarnfile1 = "genYarn.txt";
			Fiber::Yarn yarn;
			yarn.parse(configfile);	
			//yarn.setStepNum(vrtx);
			yarn.yarn_simulate();
			yarn.write_yarn(yarnfile1);
			break;
		}
		case 1: {		

			const char* yarnfile1 = "genYarn_ref.txt";
			/* This yarn is the reference yarn for shapemaching (no flyaway) */
			Fiber::Yarn yarn_ref;
			yarn_ref.parse(configfile);
			yarn_ref.yarn_simulate();
			yarn_ref.write_yarn(yarnfile1);

			std::string line;
			while (std::getline(fin00, line))
			{
				if (line == "eof") break;
				std::vector<std::string> splits = split(line, ' ');
				std::string dataset = splits[0];
				int frame0 = atoi(splits[1].c_str()) / skipFactor;
				int frame1 = atoi(splits[2].c_str()) / skipFactor + 1;
				int yarn0 = atoi(splits[3].c_str());
				int yarn1 = atoi(splits[4].c_str());
				full_pipeline(yarnfile1, configfile, vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			}

			break;
		}
		case 2: {
			std::string line;
			while (std::getline(fin00, line))
			{
				std::vector<std::string> splits = split(line, ' ');
				std::string dataset = splits[0];
				int frame0 = atoi(splits[1].c_str()) / skipFactor;
				int frame1 = atoi(splits[2].c_str()) / skipFactor + 1;
				int yarn0 = atoi(splits[3].c_str());
				int yarn1 = atoi(splits[4].c_str());
				step4_NN_output(configfile, vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			}
			break;
		}
	}

	//	std::system("pause"); //add breakpoint instead 
	return 0;
}