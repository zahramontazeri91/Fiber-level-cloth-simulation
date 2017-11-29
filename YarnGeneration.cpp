#include <fstream>
#include "Fiber.h"
#include "crossSection.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"

void fittingPlyCenter(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* plyCenterFile);
void fittingCompress(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const char* compressFile);
void fittingFiberTwisting(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* fiberTwistFile);

int main(int argc, const char **argv) {
#if 1

	if (argc != 3) {
		printf("Usage: YarnGeneration [task file] [output file]\n");
		return 1;
	}

	std::ifstream fin(argv[1]);
	if (fin.is_open()) {
		std::cout << "Using task file: \"" << argv[1] << "\"." << std::endl;

		Fiber::Yarn yarn;
		std::string command, configFILE, simulatedFILE, cntrYarnFILE,
			compressFILE, plyCenterFILE, fiberTwistFILE;

		// Fitting step 
		fin >> command >> configFILE >> simulatedFILE >> cntrYarnFILE;
		fin >> command >> plyCenterFILE >> compressFILE >> fiberTwistFILE;
		yarn.parse(configFILE.c_str());
		std::vector<yarnIntersect2D> allPlaneIntersect;
		CrossSection cs(simulatedFILE.c_str(), cntrYarnFILE.c_str(), yarn.getPlyNum(), yarn.getStepNum(), 100, allPlaneIntersect);
		fittingCompress(cs, allPlaneIntersect, compressFILE.c_str());
		fittingPlyCenter(cs, allPlaneIntersect, yarn.getYarnRadius(), plyCenterFILE.c_str());
		fittingFiberTwisting(cs, allPlaneIntersect, yarn.getYarnRadius(), fiberTwistFILE.c_str());

		// Procedural step
		yarn.yarn_simulate(plyCenterFILE.c_str(), fiberTwistFILE.c_str());
		yarn.compress_yarn(compressFILE.c_str());
		yarn.curve_yarn(cntrYarnFILE.c_str());
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


void fittingCompress(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const char* compressFile) {

	//transfer from e1-e2 to x-y plane
	std::vector<yarnIntersect2D> xy_Its;
	cs.transferLocal2XY(allPlaneIntersect, xy_Its);

	//find the fitted ellipse parameters
	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(xy_Its, ellipses);

	//regularize the ellipse params
	//std::vector<Ellipse> simple_ellipses1, simple_ellipses2, simple_ellipses3;
	//cs.regularizeEllipses(ellipses, simple_ellipses1);
	//cs.regularizeEllipses(simple_ellipses1, simple_ellipses2);

	//parameterize the ellipse params
	std::vector<Ellipse> param_ellipses;
	cs.parameterizeEllipses(ellipses, param_ellipses);

	FILE *fout;
	if (fopen_s(&fout, compressFile, "wt") == 0) {
		fprintf_s(fout, "%d \n", param_ellipses.size());
		for (int i = 0; i < param_ellipses.size(); ++i) {
			fprintf_s(fout, "%.4f %.4f %.4f \n", param_ellipses[i].longR, param_ellipses[i].shortR, param_ellipses[i].angle);
		}
		fclose(fout);
	}

}

void fittingPlyCenter(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* plyCenterFile)
{
	//fit the ellipse and find the compression param
	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect, ellipses);

	// should not simplify the ellipse params

	// Decompress simulated yarn e1-e2 space
	std::vector<yarnIntersect2D> deCompressPlaneIntersect;
	cs.deCompressYarn(allPlaneIntersect, yarn_radius, ellipses, deCompressPlaneIntersect);

	//transfer from e1-e2 to x-y plane
	std::vector<yarnIntersect2D> xy_Its;
	cs.transferLocal2XY(deCompressPlaneIntersect, xy_Its);

	//extract ply-centers helix positions
	cs.extractPlyTwist(xy_Its, plyCenterFile);

	//write plyCenters as R and theta and parameterize them
	cs.parameterizePlyCenter(plyCenterFile, plyCenterFile);
}

void fittingFiberTwisting(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* fiberTwistFile)
{
	//fit the ellipse and find the compression param
	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect, ellipses);

	//simplify the ellipse params
	std::vector<Ellipse> simple_ellipses;
	cs.parameterizeEllipses(ellipses, simple_ellipses);
	//simple_ellipses = ellipses;

	//Decompress simulated yarn e1-e2 space
	std::vector<yarnIntersect2D> deCompressPlaneIntersect;
	cs.deCompressYarn(allPlaneIntersect, yarn_radius, simple_ellipses, deCompressPlaneIntersect);

	//
	std::vector<float> fiber_theta;
	cs.fiberTwisting(deCompressPlaneIntersect, fiber_theta, fiberTwistFile);

	// no need transfer from e1-e2 to x-y plane because it's only theta

	//extract ply-centers helix positions
	//cs.extractPlyTwist(xy_Its, plyCenterFile);
}