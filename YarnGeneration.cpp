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

#if 1
		/* given a fiber.txt file, return upsampled fiber.txt */
		const int sampleRate = 10;
		std::string dataset3 = "woven/6x6";
		const int f = 6000;
		const int yarn = 12;
		for (int y = 0; y < yarn; y++) {
			//std::string tmp1 = "data/" + dataset3 + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			std::string tmp1 = "output/" + dataset3 + "/genYarn_NN_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* infile = tmp1.c_str();
			//std::string tmp2 = "data/" + dataset3 + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			std::string tmp2 = "output/" + dataset3 + "/genYarn_NN_wo_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* outfile = tmp2.c_str();
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
		return 0;
#endif

#if 1
		/* given an obj file, sample points uniformaly for specific #vertices along curve */
		/**** first generate fe and obj using "woven_generator_arbitrary.py" then use "deform.py" to extract them as usual, then run this with phase 0! ***/
		int seg_subdiv = 10;
		const int yarn0 = 0;  
		const int yarn1 = 1000;
		std::string dataset2 = "woven/arbitrary_pattern/600x400";
		for (int y = yarn0; y < yarn1; ++y) {
			std::cout << "generating for frame " << y << std::endl;

			std::string tmp1 = "input/" + dataset2 + "/centerYarn_0_" + std::to_string(y) + "_ds.txt";
			const char* curvefile = tmp1.c_str();
			std::string tmp2 = "input/" + dataset2 + "/centerYarn_sampled_0_" + std::to_string(y) + ".obj";
			const char*curvefile_new = tmp2.c_str();

			/*knitted example*/
			//const char* curvefile = "D:/sandbox/fiberSimulation/dataSets/woven/knitted/yarn/openwork_trellis_pattern.txt";
			//const char* curvefile_new = "D:/sandbox/fiberSimulation/dataSets/woven/knitted/yarn/openwork_trellis_pattern_new.obj";

			const char* normfile = "junk_norm.txt";

			HermiteCurve curve;
			curve.init(curvefile, normfile, seg_subdiv);
			double curveLength = curve.totalLength();
			double edgeLength = 0.02;
			const int pnts = int(std::ceil(curveLength / edgeLength));
			edgeLength = curveLength / static_cast<double>(pnts);
			std::ofstream fout(curvefile_new);

			for (int v = 0; v <= pnts; v++) {
				double len = v * edgeLength;
				double t = curve.arcLengthInvApprox(len);
				Eigen::Vector3d pos = curve.eval(t);
				fout << "v " << 4.f*pos[0] << " " << 4.f*pos[1] << " " << 4.f*pos[2] << std::endl;
			}

			for (int l = 1; l <= pnts; l++) {
				fout << "l " << l << " " << l + 1 << std::endl;
			}
			fout.close();
		}

		return 0;
#endif

#if 0
		/* For upsampling the stretched yarn: */ 
		//std::string dataset = "single_yarn/yarn100/stretch";
		std::string dataset = "single_yarn/yarn8/teeth/4_1.2_00110";
		const char* congif = "yarnTypes/yarn8/config_step2.txt";
		const int upsample = 3; //odd number so phase matches better
		const int vrtx = 300 * upsample;
		const int frame0 = 0;
		const int frame1 = 20000;
		const int skipfactor = 1000;
		for (int f = frame0; f < frame1 + 1; f+=skipfactor) {
			std::cout << "upsample " << f << " started .. \n";
			upsample_stretched(congif, vrtx, f, dataset, upsample);
		}
		return 0;
#endif

#if 0
		/* procedurally generate a straight yarn given the config file */
		Fiber::Yarn yarn;
		yarn.parse(argv[2]);
		yarn.simulate_ply_shuang();
		yarn.write_plys("test_ply.txt");
		const int K = yarn.getPlyNum();
		yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
		yarn.build("test_fly.txt", K);

		yarn.write_yarn("genYarn.txt");
		return 0;
#endif
	}

	if ( argc < 4 ) {
		std::cout << "Number of argument: " << argc << std::endl;
		printf("USAGE: YarnGeneration [phase1/phase2] [configFile] [datasetFile] -w window-size=5 -s upsample=2 -t isTrain=1 -x trimPercent -k skipfactor=500 -v vrtx-num -c isCompress -z stepSize_ds -rx resolution-AABB-x -ry -rz -rad radius-AABB \n");
		printf("EXAMPLE: YarnGeneration 1 yarnTypes/yarn4/config_step2.txt yarnTypes/yarn4/datasets.txt -w 5 -s 2 -t 0 -x 0.0 -k 500 -v 150 -c 1 -z 0.02 -rx 5 -ry 5 -rz 20 -rad 0.1 \n");
		return 1;
	}

	const char* configfile = argv[2];
	const char* datasetfile = argv[3];
	int window_size = 5;
	int upsample = 1;
	int isTrain = 0;
	int skipFactor = 500;
	float trimPercent = 0.0;
	int isCompress = 1;
	int vrtx_ds = 150;
	float stepSize_ds = 0.02;
	int resol_x = 1;
	int resol_y = 1;
	int resol_z = 1;
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
			std::cout << " *****      PHASE 1      ****** \n \n";
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
				int frame0 = atoi(splits[1].c_str());
				int frame1 = atoi(splits[2].c_str());
				int yarn0 = atoi(splits[3].c_str());
				int yarn1 = atoi(splits[4].c_str());
				full_pipeline(yarnfile1, configfile, vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent, upsample, stepSize);
			}
			break;
		}
		case 2: {
			std::string line;
			std::cout << " *****      PHASE 2      ****** \n \n";
			while (std::getline(fin00, line))
			{
				std::cout << " *********** \n \n";
				if (line == "eof") break;
				std::vector<std::string> splits = split(line, ' ');
				std::string dataset = splits[0];
				int frame0 = atoi(splits[1].c_str());
				int frame1 = atoi(splits[2].c_str());
				int yarn0 = atoi(splits[3].c_str());
				int yarn1 = atoi(splits[4].c_str());
				step4_NN_output(configfile, vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, isCompress, stepSize);
				step5_createVOL( skipFactor, frame0, frame1, yarn0, yarn1, dataset, resol_x, resol_y, resol_z, radius);
			}
			break;
		}
	}
	std::cout << "Code is done successfully! \n";
	//	std::system("pause"); //add breakpoint instead 
	return 0;
}