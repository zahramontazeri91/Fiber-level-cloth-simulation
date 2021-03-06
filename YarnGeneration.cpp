#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"
#include <string>
#include <cstdlib>


int main(int argc, const char **argv) {

	if (argv[0] == "-h" || argv[0] == "-help") {
		printf("USAGE: \nYarnGeneration [phase] [dataset f0 f1 y0 y1] [optional parameters] \n");
		return 1;
	}

	const int phase = atoi(argv[1]);
	assert(phase == 0 || phase == 1 || phase == 2);
	if (phase == 0 ) { //DEBUG USE

#if 0
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

//#if 0
//		// no need for this
//		/* For upsampling the stretched yarn: */ 
//		std::string dataset = "single_yarn/yarn4/stretch";
//		//std::string dataset = "single_yarn/yarn8/teeth/4_1.2_00110";
//		const char* congif = "yarnTypes/yarn4/config_step2.txt";
//		const int upsample = 3; //odd number so phase matches better
//		const int vrtx = 300 * upsample;
//		const int frame0 = 0;
//		const int frame1 = 20000;
//		const int skipfactor = 1000;
//		for (int f = frame0; f < frame1 + 1; f+=skipfactor) {
//			std::cout << "upsample " << f << " started .. \n";
//			upsample_stretched(congif, vrtx, f, dataset, upsample);
//		}
//		return 0;
//#endif

#if 0
		/* procedurally generate a straight yarn given the config file */
		Fiber::Yarn yarn;
		yarn.parse(argv[2]);
		yarn.simulate_ply_shuang();
		yarn.write_plys("test_ply.txt");
		const int K = yarn.getPlyNum();
		yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
		yarn.build("test_fly.txt", K);

		yarn.write_yarn("genYarn_30_b0_a1.txt");
		return 0;
#endif
	}

	if ( argc < 6 ) { //first six arguments are necessary: [phase dataset f0 f1 y0 y1] -optional params
		std::cout << "Number of argument: " << argc << std::endl;
		printf("USAGE: \nYarnGeneration [phase1/phase2] [dataset f0 f1 y0 y1] -w window-size=5 -s upsample=2 -t isTrain=1 -x trimPercent -v vrtx-num -c isCompress -z stepSize_ds -v hasVol -rx resolution-AABB-x -ry -rz -rad radius-AABB \n");
		printf("EXAMPLE: \nYarnGeneration 1 yarn8/pattern/spacing0.5x/10 7000 14000 0 1 -k 500 -v 150 -t 0 -s 2 -w 5 -s 2 -t 0 -x 0.0 -k 500 -v 150 -c 1 -scale 0.25 -z 0.02 -vol 1 -rx 5 -ry 5 -rz 20 -rad 0.1 \n");
		return 1;
	}

	std::string dataset_str = argv[2];
	std::string f0_str = argv[3];
	std::string f1_str = argv[4];
	std::string y0_str = argv[5];
	std::string y1_str = argv[6];

	const char* dataset = dataset_str.c_str();
	int indx = dataset_str.find('/');
	std::string yarnType = dataset_str.substr(0, indx);
	std::string config_str = "yarnTypes/" + yarnType + "/config_step2.txt";
	const char* configfile = config_str.c_str();

	int frame0 = atoi(f0_str.c_str());
	int frame1 = atoi(f1_str.c_str());
	int yarn0 = atoi(y0_str.c_str());
	int yarn1 = atoi(y1_str.c_str());

	// Default parameters:
	int window_size = 5;
	int upsample = 2;
	int upsampleMore = 1;
	int isTrain = 0;
	float trimPercent = 0.1f;
	int isCompress = 1;
	int vrtx_ds = 150;
	float stepSize_ds = 0.02f;
	int hasVol = 0;
	int resol_x = 1;
	int resol_y = 1;
	int resol_z = 1;
	float radius = 0.1f;
	float scaleSim = 0.25f;
	int skipFactor = 500;

	for (int i = 5; i < argc-1; ++i) { //first six args are already read
		std::string arg = argv[i];
		std::stringstream arg1(argv[i+1]);
		if (arg == "-k")
			arg1 >> skipFactor;
		if (arg == "-w")
			arg1 >> window_size;
		if (arg == "-s")
			arg1 >> upsample;
		if (arg == "-ss")
			arg1 >> upsampleMore;
		if (arg == "-t")
			arg1 >> isTrain;
		if (arg == "-x")
			arg1 >> trimPercent;
		if (arg == "-v")
			arg1 >> vrtx_ds;
		if (arg == "-c")
			arg1 >> isCompress;
		if (arg == "-scale")
			arg1 >> scaleSim;
		if (arg == "-z")
			arg1 >> stepSize_ds;
		if (arg == "-vol")
			arg1 >> hasVol;
		if (arg == "-rx")
			arg1 >> resol_x;
		if (arg == "-ry")
			arg1 >> resol_y;
		if (arg == "-rz")
			arg1 >> resol_z;
		if (arg == "-rad")
			arg1 >> radius;
	}
	const float stepSize = stepSize_ds / float(upsample);
	const int vrtx = vrtx_ds * upsample;

	//if the data is for training, both ends must have some trimPercent
	if (isTrain)
		assert(trimPercent > 0 && "For training, both ends must be trimmed, choose -t larger than 0");

	std::ifstream fin0(configfile);
	//std::ifstream fin00(datasetfile);
	if (phase) {
		assert(fin0.is_open() && "config file wasn't found!\n");
	}

	switch (phase) {
		case 1: {		
			std::cout << " *****      PHASE 1      ****** \n \n";
			generateNNinput(configfile, vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, scaleSim, window_size, trimPercent, upsample, yarnType);
			break;
		}

		case 2: {
			std::cout << " *****      PHASE 2      ****** \n \n";

			/** Given NN output, deform the procedural yarn **/
			step5_applyNNoutput(configfile, vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, isCompress, stepSize);
				
			/** optional step of generating bounding box needed for mitsuba-ct **/
			if (hasVol)
				step6_createVOL( skipFactor, frame0, frame1, yarn0, yarn1, dataset, resol_x, resol_y, resol_z, radius);

			/** optional step for upsampling stretched yarns **/
			if (upsampleMore != 1)
				step7_upsample(skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, isCompress, upsampleMore);
			break;
		}
	}
 	std::cout << "Code is done successfully! \n";
	//	std::system("pause"); //add breakpoint instead 	return 0;
}