#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"

int main(int argc, const char **argv) {
	int phase = 1;
	//phase 0: test
	//phase 1: fitting
	//phase 2: training
	//phase 3: parameterization 
	//phase 4: generation

	switch (phase) {
		case 1: {
			std::cout << "*** Fitting phase ***\n";
			const char* configfile = "config.txt";

			const char* yarnfile1 = "data/frame0000000.txt";
			//const char* yarnfile1 = "data/genYarn_frame1_shuang.txt";
			//const char* yarnfile1 = "data/frame00001_scaled.txt";
			//const char* yarnfile1 = "data/genYarn_frame1.txt"; 

			const char* yarnfile2 = "data/frame0011100.txt";
			//const char* yarnfile2 = "data/frame00001_scaled.txt";
			//const char* yarnfile2 = "data/frame00029_scaled.txt";
			//const char* yarnfile2 = "data/genYarn_frame1_compressed_R.txt";
			//const char* norm2 = "data/genYarn_frame29_norms.txt";

			const char* compressfile = "compressParams.txt";

			std::ifstream fin1(configfile);
			std::ifstream fin2(yarnfile1);
			std::ifstream fin3(yarnfile2);
			assert(fin1.is_open() && "config file wasn't found!\n");
			assert(fin2.is_open() && "compressed-yarn file wasn't found!\n");
			assert(fin3.is_open() && "reference-yarn file wasn't found!\n");

			extractCompressParams(configfile, yarnfile1, yarnfile2, compressfile);
			break;
		}
		case 2: {
			std::cout << "*** Training phase ***\n";
			break;
		}
		case 3: {
			std::cout << "*** Parameterization phase ***\n";
			break;
		}
		case 4: {
			std::cout << "*** Generation phase *** \n";

			Fiber::Yarn yarn;
			const char* configfile = "config.txt";
			const char* compressfile = "compressParams.txt";
			const char* curvefile = "centerYarn_compress.txt";

			std::ifstream fin1(configfile);
			std::ifstream fin2(compressfile);
			std::ifstream fin3(curvefile);
			assert(fin1.is_open() && "config file wasn't found!\n");
			assert(fin2.is_open() && "compressfile file wasn't found!\n");
			assert(fin3.is_open() && "curvefile file wasn't found!\n");

			yarn.parse(configfile);

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