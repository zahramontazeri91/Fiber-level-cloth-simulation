#include <fstream>
#include "Fiber.h"
#include "crossSection.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"

void fittingPlyCenter(CrossSection & cs, const char* plyCenterFile);
void fittingCompress(CrossSection & cs, const char* compressFile, const char* normsFile);
void fittingFiberTwisting(CrossSection & cs, const char* fiberTwistFile);

int main(int argc, const char **argv) {
#if 1

    if ( argc != 3 ) {
        printf("Usage: YarnGeneration [task file] [output file]\n");
		return 1;
    }

    std::ifstream fin(argv[1]);
	if (fin.is_open()) {
		std::cout << "Using task file: \"" << argv[1] << "\"." << std::endl;

		Fiber::Yarn yarn;
		std::string command, configFILE, simulatedFILE,
			cntrYarnFILE, curveNorms, compressFILE, plyCenterFILE;

		fin >> configFILE;
		yarn.parse(configFILE.c_str());

		// Fitting step 
		fin >> command >> simulatedFILE >> cntrYarnFILE >> plyCenterFILE >> compressFILE >> curveNorms;
		if (command == "FITTING") {
			CrossSection cs(simulatedFILE.c_str(), cntrYarnFILE.c_str(), yarn.getPlyNum(), yarn.getStepNum(), 100);
			fittingCompress(cs, compressFILE.c_str(), curveNorms.c_str());
			fittingPlyCenter(cs, plyCenterFILE.c_str());
			fittingFiberTwisting(cs, "fiberTwist.txt");
		}
		// Procedural step
		fin >> command >> plyCenterFILE;
		if (command == "SIMULATE")
			yarn.yarn_simulate(plyCenterFILE.c_str(), "fiberTwist.txt");
		//uncomment if generate yarn without given ply-center
			//yarn.yarn_simulate();

		fin >> command >> compressFILE;
		if (command == "COMPRESS")
			yarn.compress_yarn(compressFILE.c_str());

		fin >> command >> curveNorms;
		if (command == "CURVE")
			yarn.curve_yarn(cntrYarnFILE.c_str(), curveNorms.c_str());

		yarn.write_yarn(argv[2]);
	}
	else
		std::cout << "File wasn't found! \n";


#else
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
	ply_centers_test();

	//extractNormals();
#endif

	//std::system("pause"); //add breakpoint instead

	return 0;
}


void fittingCompress(CrossSection & cs, const char* compressFile, const char* normsFile) {
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);

	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);

	//transfer from e1-e2 to x-y plane
	std::vector<yarnIntersect2D> xy_Its;
	cs.transferLocal2XY(allPlaneIntersect, xy_Its);

	//find the fitted ellipse parameters
	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(xy_Its, ellipses);

	//simplify the ellipse params
	std::vector<Ellipse> simple_ellipses;
	//cs.parameterizeEllipses(ellipses, simple_ellipses);
	simple_ellipses = ellipses;

	FILE *fout;
	if (fopen_s(&fout, compressFile, "wt") == 0) {
		fprintf_s(fout, "%d \n", simple_ellipses.size());
		for (int i = 0; i < simple_ellipses.size(); ++i) {
			fprintf_s(fout, "%.4f %.4f %.4f \n", simple_ellipses[i].longR, simple_ellipses[i].shortR, simple_ellipses[i].angle);
		}
		fclose(fout);
	}

	//extract spline normals
	//std::vector<vec3f> normals;
	//cs.extractNormals(ellipses, normals, normsFile);
}

void fittingPlyCenter(CrossSection & cs, const char* plyCenterFile )
{
	//find 3D intersections with plane
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);

	//project 3D intersections in e1-e2 plane
	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);

	//fit the ellipse and find the compression param
	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect, ellipses);

	//simplify the ellipse params
	std::vector<Ellipse> simple_ellipses;
	cs.parameterizeEllipses(ellipses, simple_ellipses);
	//simple_ellipses = ellipses;

	// Decompress simulated yarn e1-e2 space
	std::vector<yarnIntersect2D> deCompressPlaneIntersect;
	cs.deCompressYarn(allPlaneIntersect, simple_ellipses, deCompressPlaneIntersect);

	//transfer from e1-e2 to x-y plane
	std::vector<yarnIntersect2D> xy_Its;
	cs.transferLocal2XY(deCompressPlaneIntersect, xy_Its);

	//extract ply-centers helix positions
	cs.extractPlyTwist(xy_Its, plyCenterFile);
}

void fittingFiberTwisting(CrossSection & cs, const char* fiberTwistFile)
{
	//find 3D intersections with plane
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);

	//project 3D intersections in e1-e2 plane
	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);

	//fit the ellipse and find the compression param
	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect, ellipses);

	//simplify the ellipse params
	std::vector<Ellipse> simple_ellipses;
	//cs.parameterizeEllipses(ellipses, simple_ellipses);
	simple_ellipses = ellipses;

	//Decompress simulated yarn e1-e2 space
	std::vector<yarnIntersect2D> deCompressPlaneIntersect;
	cs.deCompressYarn(allPlaneIntersect, simple_ellipses, deCompressPlaneIntersect);

	//
	std::vector<float> fiber_theta;
	cs.fiberTwisting(deCompressPlaneIntersect, fiber_theta, fiberTwistFile);

	// no need transfer from e1-e2 to x-y plane because it's only theta

	//extract ply-centers helix positions
	//cs.extractPlyTwist(xy_Its, plyCenterFile);
}