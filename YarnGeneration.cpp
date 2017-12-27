#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"
#include <string>

int main(int argc, const char **argv) {

	const char* configfile = "config.txt";
	std::ifstream fin1(configfile);
	assert(fin1.is_open() && "config file wasn't found!\n");
	Fiber::Yarn yarn;
	yarn.parse(configfile);
	
	int frame0 = 35;
	int frame1 = 36;
	std::string dataset = "1220";
	//std::string tmp0 = "data/" + dataset + "/simul_frame_0.txt";
	//const char* yarnfile1 = tmp0.c_str();
	const char* yarnfile1 = "genYarn_ref.txt";

	//phase 0: test
	//phase 1: fitting
	//phase 2: training
	//phase 3: parameterization 
	//phase 4: curved-yarn generation

	int phase = 1;


	switch (phase) {
		case 1: {
			std::cout << "*** Fitting phase ***\n";

			for (int i = frame0; i < frame1; i++) {		

				std::string tmp1 = "data/" + dataset + "/simul_frame_" + std::to_string(i * 5) + ".txt";
				const char* yarnfile2 = tmp1.c_str();
				std::string tmp2 = "input/" + dataset + "/matrix_R_" + std::to_string(i * 5) + ".txt";
				const char* compress_R = tmp2.c_str();
				std::string tmp3 = "input/" + dataset + "/matrix_S_" + std::to_string(i * 5) + ".txt";
				const char* compress_S = tmp3.c_str();
				std::string tmp4 = "input/" + dataset + "/centerYarn_" + std::to_string(i * 5) + ".txt";
				const char* curvefile = tmp4.c_str();
				std::string tmp5 = "input/" + dataset + "/normYarn_" + std::to_string(i * 5) + ".txt";
				const char* normfile = tmp5.c_str();

				std::ifstream fin1(yarnfile1);
				std::ifstream fin2(yarnfile2);
				assert(fin1.is_open() && "reference-yarn file wasn't found!\n");
				assert(fin2.is_open() && "compressed-yarn file wasn't found!\n");

				const int vrtx_num = yarn.getStepNum();
				extractCompress_seg(yarnfile1, yarnfile2, compress_R, compress_S, 
					curvefile, normfile, yarn.getPlyNum(), vrtx_num);

				//Fiber::Yarn::Compress compress;
				//Fiber::Yarn::CenterLine curve;
				//const float trimPercent = 0.33;
				//constFitting_compParam(ellipses, theta_R, trimPercent, compress);
				//sinFitting_curve(curvefile, trimPercent, curve);

				/**************************************************/
				std::string tmp6 = "output/" + dataset + "/genYarn_" + std::to_string(i * 5) + ".txt";
				const char* outfile = tmp6.c_str();
				// Procedural step
				yarn.yarn_simulate();
				yarn.compress_yarn(compress_R, compress_S);
				yarn.curve_yarn(curvefile, normfile);
				yarn.write_yarn(outfile);

				//std::string tmp7 = "output/" + dataset + "/genYarn_wo_" + std::to_string(i * 5) + ".txt";
				//const char* outfile_wo = tmp7.c_str();
				//yarn.yarn_simulate();
				////yarn.compress_yarn(compress_R, compress_S);
				//yarn.curve_yarn(curvefile, normfile);
				//yarn.write_yarn(outfile_wo);
			}
			break;
		}
		case 2: {
			std::cout << "*** Training phase ***\n";

			for (int i = frame0; i < frame1; i++) {

				std::string tmp5 = "input/" + dataset + "/matrix_R_" + std::to_string(i * 5) + ".txt";
				const char* compress_R = tmp5.c_str();
				std::string tmp6 = "input/" + dataset + "/NN/testY_NN_full" + std::to_string(i * 5) + ".txt";
				const char* compress_S = tmp6.c_str();
				std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(i * 5) + ".txt";
				const char* curvefile = tmp7.c_str();
				std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(i * 5) + ".txt";
				const char* normfile = tmp8.c_str();


				std::string tmp3 = "output/" + dataset + "/genYarn_NN_" + std::to_string(i * 5) + ".txt";
				const char* outfile = tmp3.c_str();
				// Procedural step
				yarn.yarn_simulate();
				yarn.compress_yarn(compress_R, compress_S);
				yarn.curve_yarn(curvefile, normfile);
				yarn.write_yarn(outfile);
			}

			break;
		}
		case 3: {
			//std::cout << "*** Parameterization and merging phase ***\n";
			//const char* compressFile = "compressParam_yarn.txt";
			//const char* curveFile = "centerline_yarn.txt";
			//std::vector<Fiber::Yarn::Compress> compress_segs;
			//const int seg_vrtx = 150;
			//const int seg_num = 3;
			//const int yarn_vrtx = yarn.getStepNum();

			//compress_segs.resize(seg_num);
			//for (int i = 0; i < seg_num; ++i) {
			//	compress_segs[i].ellipse_long = 0.877276;
			//	compress_segs[i].ellipse_short = 0.452084;
			//	compress_segs[i].ellipse_theta = 2.43881;
			//	compress_segs[i].rotation = 5.23583;
			//}
			//appendCompress_yarn(compress_segs, seg_vrtx, yarn_vrtx, compressFile );

			//std::vector<Fiber::Yarn::CenterLine> centerlines;

			//centerlines.resize(seg_num);
			//for (int i = 0; i < seg_num; ++i) {
			//	centerlines[i].a = 0.0091413;
			//	centerlines[i].b = 0.12083;
			//	centerlines[i].c = 0.00551934;
			//}

			//appendCenter_yarn(centerlines, seg_vrtx, yarn_vrtx, curveFile);
			
			break;
		}
		case 4: {
			std::cout << "*** Generation phase *** \n";
			//Generate a yarn mapped to a given curve
			const char* configFile = argv[1];
			std::ifstream fin1(configFile);
			assert(fin1.is_open() && "config file wasn't found!\n");

			Fiber::Yarn yarn0;
			yarn0.parse(configFile);

			const char* curvefile = argv[2];
			//const char* normfile = "normYarn.txt";
			std::ifstream fin2(curvefile);
			assert(fin2.is_open() && "curve file wasn't found!\n");

			// Procedural step
			yarn0.yarn_simulate();
			yarn0.curve_yarn(curvefile);
			yarn0.write_yarn(argv[3]);
			break;
		}
		case 0: {
			std::cout << "*** Testing ***\n";
			//hermiteTest1();
			//hermiteTest2();
			//hermiteTest3();

			//linePlaneIntersection_test();
			//yarnPlaneIntersection_test();
			//bildPlanes_test();

			//allPlanesIntersections_test(); 
			//project2Plane_test();
			//write_PlanesIntersections2D_test();
			//getOrientation_test();

			//extractCompressParam_test();

			//compress_yarn_test();
			//ply_centers_test();

			//extractNormals();

			//shapeMatch_test();
			//yarnShapeMatch_test();

			//writeNormals();

			//render1fiber();
			break;
		}
	}

	//	std::system("pause"); //add breakpoint instead
 	return 0;
}