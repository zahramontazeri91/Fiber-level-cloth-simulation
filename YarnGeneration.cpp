#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"

int main(int argc, const char **argv) {
	const char* configfile = "config.txt";
	std::ifstream fin1(configfile);
	assert(fin1.is_open() && "config file wasn't found!\n");
	Fiber::Yarn yarn;
	yarn.parse(configfile);

	int phase = 3;
	//phase 0: test
	//phase 1: fitting
	//phase 2: training
	//phase 3: parameterization 
	//phase 4: generation

	switch (phase) {
		case 1: {
			std::cout << "*** Fitting phase ***\n";

			//const char* yarnfile1 = "data/frame0000000.txt";
			const char* yarnfile1 = "data/genYarn_frame1_shuang.txt";
			//const char* yarnfile1 = "data/frame00001_scaled.txt";
			//const char* yarnfile1 = "data/genYarn_frame1.txt"; 

			//const char* yarnfile2 = "data/frame0011100.txt";
			//const char* yarnfile2 = "data/frame00001_scaled.txt";
			const char* yarnfile2 = "data/frame00029_scaled.txt";
			//const char* yarnfile2 = "data/genYarn_frame1_compressed_R.txt";
			//const char* norm2 = "data/genYarn_frame29_norms.txt";

			const char* compressfile = "compressParams_seg.txt";

			std::ifstream fin1(yarnfile1);
			std::ifstream fin2(yarnfile2);
			assert(fin1.is_open() && "compressed-yarn file wasn't found!\n");
			assert(fin2.is_open() && "reference-yarn file wasn't found!\n");

			//const int vrtx_num = yarn.getStepNum();
			const int vrtx_num = 150;
			extractCompress_seg(yarnfile1, yarnfile2, compressfile, yarn.getPlyNum(), vrtx_num);
			break;
		}
		case 2: {
			std::cout << "*** Training phase ***\n";
			break;
		}
		case 3: {
			std::cout << "*** Parameterization phase ***\n";
			const char* compressFile = "compressParam_yarn.txt";
			const char* curveFile = "centerline_yarn.txt";
			std::vector<Fiber::Yarn::Compress> compress_segs;
			const int seg_vrtx = 150;
			const int seg_num = 3;
			const int yarn_vrtx = yarn.getStepNum();
			
			compress_segs.resize(seg_num);
			for (int i = 0; i < seg_num; ++i) {
				compress_segs[i].ellipse_long = 1.7;
				compress_segs[i].ellipse_short = 0.3;
				compress_segs[i].ellipse_theta = 0.0;
				compress_segs[i].rotation = 0.0;
			}
			appendCompress_yarn(compress_segs, seg_vrtx, yarn_vrtx, compressFile );

			std::vector<Fiber::Yarn::CenterLine> centerlines;

			centerlines.resize(seg_num);
			for (int i = 0; i < seg_num; ++i) {
				centerlines[i].a = 1.7;
				centerlines[i].b = 0.3;
				centerlines[i].c = 0.0;
				centerlines[i].d = 1.0;
			}

			appendCenter_yarn(centerlines, seg_vrtx, yarn_vrtx, curveFile);
			
			break;
		}
		case 4: {
			std::cout << "*** Generation phase *** \n";

			const char* compressfile = "compressParams_seg.txt";
			const char* curvefile = "centerYarn_compress.txt";

			std::ifstream fin1(compressfile);
			std::ifstream fin2(curvefile);
			assert(fin1.is_open() && "compressfile file wasn't found!\n");
			assert(fin2.is_open() && "curvefile file wasn't found!\n");

			// Procedural step
			yarn.yarn_simulate();
			yarn.compress_yarn(compressfile);
			yarn.curve_yarn(compressfile, curvefile);
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