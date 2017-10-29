#include <fstream>
#include "Fiber.h"
#include "crossSection.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"

void fittingPlyCenter(CrossSection & cs, const char* compressFile, const char* plyCenterFile);
void fittingCompress(CrossSection & cs, const char* compressFile);
int main(int argc, const char **argv) {
#if 1

    if ( argc != 3 ) {
        printf("Usage: YarnGeneration [task file] [output file]\n");
		std::system("pause");
		return 1;
    }

    std::ifstream fin(argv[1]);
	if (fin.is_open()) {
		std::cout << "Using task file: \"" << argv[1] << "\"." << std::endl;

		Fiber::Yarn yarn;
		std::string command, configFILE, simulatedFILE, 
			curveFILE, compressFILE, plyCenterFILE;

		fin >> configFILE;
		yarn.parse(configFILE.c_str());

		// Fitting step 
		fin >> command >> simulatedFILE >> curveFILE >> compressFILE >> plyCenterFILE;
		if (command == "FITTING") {
			CrossSection cs(simulatedFILE.c_str(), curveFILE.c_str(), yarn.getPlyNum(), yarn.getStepNum(), 100);
			fittingCompress(cs, compressFILE.c_str() );
			fittingPlyCenter(cs, compressFILE.c_str(), plyCenterFILE.c_str() );
		}
		// Procedural step
		yarn.yarn_simulate(plyCenterFILE.c_str());
		//yarn.yarn_simulate();

		fin >> command >> compressFILE;
		if (command == "COMPRESS")
			yarn.compress_yarn(compressFILE.c_str());

		fin >> command >> curveFILE;
		if (command == "CURVE")
			yarn.curve_yarn(curveFILE.c_str());

		yarn.write_yarn(argv[2]);
	}
	else
		std::cout << "File wasn't found! \n";


#else
    //hermiteTest1();
    //hermiteTest2();

	//linePlaneIntersection_test();
	//yarnPlaneIntersection_test();
	//bildPlanes_test();
	//allPlanesIntersections_test(); 
	//project2Plane_test();
	//write_PlanesIntersections2D_test();
	//getOrientation_test();
	//extractCompressParam_test();

	//compress_yarn_test();
	ply_centers_test();
#endif

	//std::system("pause"); //add breakpoint instead
	return 0;
}

void fittingCompress(CrossSection & cs, const char* compressFile) {
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);

	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);

	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect, ellipses, compressFile);
}

void fittingPlyCenter(CrossSection & cs, const char* compressFile, const char* plyCenterFile )
{
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);
	//TODO: move cs to main() and only pass it to these two functions
	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect, ellipses, compressFile);
	// Decompress simulated yarn
	std::vector<yarnIntersect2D> deCompressPlaneIntersect;
	cs.deCompressYarn(allPlaneIntersect, ellipses, deCompressPlaneIntersect);

	//extract ply-centers helix positions
	cs.extractPlyTwist(deCompressPlaneIntersect, plyCenterFile);
}