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

	int phase = 1;
	//phase 0: test
	//phase 1: fitting
	//phase 2: training
	//phase 3: parameterization 
	//phase 4: generation

	switch (phase) {
		case 1: {
			std::cout << "*** Fitting phase ***\n";

			const char* yarnfile1 = "data/1218_frame_0006000fiber_00.txt";
			std::string s = "";
			for (int i = 7; i < 8; i++) {
				if (i < 10)
					s = "frame_000" + std::to_string(i) + "000";
				else 
					s = "frame_00" + std::to_string(i) + "000";

				//const char* yarnfile2 = ("data/1218_" + s + "fiber_00.txt").c_str();
				//const char* compressfile = ("NN/param_" + s + ".txt").c_str();
				//const char* curvefile = ("centerYarn_" + s + ".txt").c_str();

				const char* yarnfile2 = "data/1218_frame_0015000fiber_00.txt";
				const char* compressfile_vrtx = "NN/param_frame_0015000.txt";
				const char* curvefile = "centerYarn_frame_0015000.txt";
				const char* compressfile_seg = "compressParams.txt";

				std::ifstream fin1(yarnfile1);
				std::ifstream fin2(yarnfile2);
				assert(fin1.is_open() && "reference-yarn file wasn't found!\n");
				assert(fin2.is_open() && "compressed-yarn file wasn't found!\n");


				//const int vrtx_num = yarn.getStepNum();
				const int vrtx_num = 150;
				std::vector<Ellipse> ellipses;
				std::vector<float> theta_R;
				extractCompress_seg(yarnfile1, yarnfile2, compressfile_vrtx, compressfile_seg, curvefile, yarn.getPlyNum(), vrtx_num, ellipses, theta_R);

				//Fiber::Yarn::Compress compress;
				//Fiber::Yarn::CenterLine curve;
				//const float trimPercent = 0.33;
				//constFitting_compParam(ellipses, theta_R, trimPercent, compress);
				//sinFitting_curve(curvefile, trimPercent, curve);
			}
			break;
		}
		case 2: {
			std::cout << "*** Training phase ***\n";
			break;
		}
		case 3: {
			std::cout << "*** Parameterization and merging phase ***\n";
			const char* compressFile = "compressParam_yarn.txt";
			const char* curveFile = "centerline_yarn.txt";
			std::vector<Fiber::Yarn::Compress> compress_segs;
			const int seg_vrtx = 150;
			const int seg_num = 3;
			const int yarn_vrtx = yarn.getStepNum();

			compress_segs.resize(seg_num);
			for (int i = 0; i < seg_num; ++i) {
				compress_segs[i].ellipse_long = 0.877276;
				compress_segs[i].ellipse_short = 0.452084;
				compress_segs[i].ellipse_theta = 2.43881;
				compress_segs[i].rotation = 5.23583;
			}
			appendCompress_yarn(compress_segs, seg_vrtx, yarn_vrtx, compressFile );

			std::vector<Fiber::Yarn::CenterLine> centerlines;

			centerlines.resize(seg_num);
			for (int i = 0; i < seg_num; ++i) {
				centerlines[i].a = 0.0091413;
				centerlines[i].b = 0.12083;
				centerlines[i].c = 0.00551934;
			}

			appendCenter_yarn(centerlines, seg_vrtx, yarn_vrtx, curveFile);
			
			break;
		}
		case 4: {
			std::cout << "*** Generation phase *** \n";

			const char* compressfile = "compressParam_yarn.txt";
			const char* curvefile = "frame00029_avg.txt";
			const char* normfile = "frame00029_norms.txt";

			std::ifstream fin1(compressfile);
			std::ifstream fin2(curvefile);
			std::ifstream fin3(normfile);
			assert(fin1.is_open() && "compressfile file wasn't found!\n");
			assert(fin2.is_open() && "curvefile file wasn't found!\n");
			assert(fin3.is_open() && "normfile file wasn't found!\n");

			// Procedural step
			yarn.yarn_simulate();
			//yarn.compress_yarn(compressfile);
			yarn.curve_yarn(curvefile, normfile);
			yarn.write_yarn(argv[2]);
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

			render1fiber();
			break;
		}
	}

	//	std::system("pause"); //add breakpoint instead
 	return 0;
}