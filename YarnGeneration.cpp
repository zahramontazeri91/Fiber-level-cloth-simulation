#include <fstream>
#include "Fiber.h"
#include "crossSection.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"

void fittingPlyCenter(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* plyCenterFile);
void fittingCompress(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const char* compressFile);
void fittingFiberTwisting(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const float yarn_radius, const char* fiberTwistFile);
void plotIntersections(const std::vector<yarnIntersect2D> &its, const char* filename);
void plotIntersectionsDeformed(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Ellipse> &ellipses, const std::vector<float> &all_theta_R, const char* filename);

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
			normFILE, compressFILE, plyCenterFILE, fiberTwistFILE;

		// Fitting step 
		fin >> command >> configFILE >> simulatedFILE >> cntrYarnFILE >> normFILE;
		fin >> command >> plyCenterFILE >> compressFILE >> fiberTwistFILE;
		yarn.parse(configFILE.c_str());
		std::vector<yarnIntersect2D> allPlaneIntersect;
		CrossSection cs(simulatedFILE.c_str(), cntrYarnFILE.c_str(), normFILE.c_str(), yarn.getPlyNum(), yarn.getStepNum(), 100, allPlaneIntersect);
		fittingCompress(cs, allPlaneIntersect, compressFILE.c_str());
		//fittingPlyCenter(cs, allPlaneIntersect, yarn.getYarnRadius(), plyCenterFILE.c_str());
		//fittingFiberTwisting(cs, allPlaneIntersect, yarn.getYarnRadius(), fiberTwistFILE.c_str());

		// Procedural step
		yarn.yarn_simulate(plyCenterFILE.c_str(), fiberTwistFILE.c_str());
		Fiber::Yarn::compress params;
		params.ellipse_long = 2.0;
		params.ellipse_short = 0.5;
		params.ellipse_theta = 0.0;
		params.rotation = pi;

		yarn.compress_yarn(params);
		//yarn.compress_yarn(compressFILE.c_str());

		//map external yarn to the curve
		//yarn.build("genYarn_frame1_compressed.txt", 2);
		
		yarn.curve_yarn(cntrYarnFILE.c_str(), normFILE.c_str());
		yarn.write_yarn(argv[2]);
	}
	else
		std::cout << "File wasn't found! \n";

#elif 0
	//if (argc != 4) {
	//	printf("Usage: YarnGeneration [config file] [reference-yarn file] [compressed-yarn file]\n");
	//	return 1;
	//}
	//const char* configfile = argv[1];
	const char* configfile = "config.txt";

	//const char* yarnfile1 = argv[2];
	const char* yarnfile1 = "frame0000000.txt";
	//const char* yarnfile1 = "genYarn_frame1_shuang.txt";
	//const char* yarnfile1 = "genYarn_frame1.txt"; 

	//const char* yarnfile2 = argv[3];
	const char* yarnfile2 = "frame0009900.txt";
	//const char* yarnfile2 = "frame00001_scaled.txt";
	//const char* yarnfile2 = "frame00029_scaled.txt";
	//const char* yarnfile2 = "genYarn_frame1_compressed_R.txt";
	//const char* norm2 = "genYarn_frame29_norms.txt";

	std::ifstream fin1(configfile);
	std::ifstream fin2(yarnfile1);
	std::ifstream fin3(yarnfile2);
	assert(fin1.is_open() && "config file wasn't found!\n");
	assert(fin2.is_open() && "compressed-yarn file wasn't found!\n");
	assert(fin3.is_open() && "reference-yarn file wasn't found!\n");

	Fiber::Yarn yarn;

	yarn.parse(configfile);
	//const int n = yarn.getStepNum();
	const int n = 150;
	const int ply_num = yarn.getPlyNum();

	Fiber::Yarn yarn_tmp;
	const char* centerYarn1 = "centerYarn_ref.txt";
	yarn_tmp.yarnCenter(yarnfile1, centerYarn1);
	std::vector<yarnIntersect2D> pnts_ref;
	CrossSection cs2(yarnfile1, centerYarn1, ply_num, n, 100, pnts_ref);

	const char* centerYarn2 = "centerYarn_compress.txt";
	yarn_tmp.yarnCenter(yarnfile2, centerYarn2);
	std::vector<yarnIntersect2D> pnts_trans;
	CrossSection cs(yarnfile2, centerYarn2, ply_num, n, 100, pnts_trans);

	//0. yarnShapeMatches:
	std::cout << "\n";
	std::cout << "Step 1: Shape matching between reference and compressed yarn.\n";
	std::vector<Ellipse> ellipses;
	std::vector<float> all_theta_R;
	cs.yarnShapeMatches(pnts_trans, pnts_ref, ellipses, all_theta_R);
	std::vector<bool> isValid(n, true);

	/*use PCA */
	//0. extract the ellipses 
	//std::vector<float> all_theta_R (n, 0.f);
	//std::vector<Ellipse> ellipses1(n);
	//std::vector<bool> isValid(n, true);
	//cs.fitEllipses(pnts_ref, ellipses1, isValid);
	//std::vector<Ellipse> ellipses2(n);
	//cs.fitEllipses(pnts_trans, ellipses2, isValid);
	//std::vector<Ellipse> ellipses(n);
	//for (int i = 0; i < n; ++i) {
	//	Ellipse ell;
	//	ell.longR = ellipses2[i].longR / ellipses1[i].longR;
	//	ell.shortR = ellipses2[i].shortR / ellipses1[i].shortR;
	//	ell.angle = ellipses2[i].angle;
	//	ellipses[i] = ell;
	//}

	//1. precompute R1\R2\theta and isValid
	std::cout << "Step 2: precompute four possibilites for each cross-section.\n";
	Eigen::MatrixXf R1(n, 4);
	Eigen::MatrixXf R2(n, 4);
	Eigen::MatrixXf theta(n, 4);
	cs.preComputeEllipses(ellipses, R1, R2, theta);

	//2. precompute cost function 
	std::cout << "Step 3: computing the cost function...\n";
	std::vector<Eigen::Matrix4f> cost;
	cs.costFunction(R1, R2, theta, isValid, cost);

	//3. dynamic programming
	std::cout << "Step 4: Use Dynamic Programming to find the optimum.\n";
	Eigen::MatrixXf totalCost(n, 4);
	Eigen::MatrixXf preConfig(n, 4);
	cs.dynamicProgramming(isValid, cost, totalCost, preConfig);

	//4.retreive solution for valid cross-sections
	std::cout << "Step 5: finalize optimum solution.\n";
	std::vector<Ellipse> validEllipses(n);
	cs.retreiveSol(R1, R2, theta, totalCost, preConfig, isValid, validEllipses);

	//5. regularizing the parameters
	//std::cout << "Step 6: regularizing the parameters.\n\n";
	//std::vector<Ellipse> simple_ellipses;
	//std::vector<float> simple_theta_R;
	//cs.regularizeEllipses(validEllipses, all_theta_R, simple_ellipses, simple_theta_R, 40);
	//cs.regularizeEllipses(simple_ellipses, simple_theta_R, validEllipses, all_theta_R, 40);

	std::cout << "\n";

	FILE *fout;
	if (fopen_s(&fout, "compressParams.txt", "wt") == 0) {
		fprintf_s(fout, "%d \n", ellipses.size());
		for (int i = 0; i < ellipses.size(); ++i) {
			//fprintf_s(fout, "%.6f %.6f %.6f ", validEllipses[i].longR, validEllipses[i].shortR, validEllipses[i].angle);
			//fprintf_s(fout, "%.6f %.6f %.6f \n", all_theta_R[i], all_T[i](0, 0), all_T[i](1, 0));
			fprintf_s(fout, "%.6f %.6f %.6f %.6f \n", validEllipses[i].longR, validEllipses[i].shortR, validEllipses[i].angle, all_theta_R[i]);
		}
		fclose(fout);
	}
	std::cout << "Compression parameters are successfully written to the file!\n";

	//visualization

	const char* refFile = "../data/allCrossSection2D_ref.txt";
	const char* deformedRefFile = "../data/allCrossSection2D_deformedRef.txt";
	const char* deformed = "../data/allCrossSection2D_deformed.txt";
	plotIntersections(pnts_ref, refFile);
	std::vector<yarnIntersect2D> ref_deformed;
	plotIntersectionsDeformed(pnts_ref, ref_deformed, ellipses, all_theta_R, deformedRefFile);
	plotIntersections(pnts_trans, deformed);
#elif 0

	Fiber::Yarn yarn;
	const char* configfile = "config.txt";
	yarn.parse(configfile);
	const int n = yarn.getStepNum();
	const int ply_num = yarn.getPlyNum();

	for (int i = 10; i < 30; ++i) {
		std::string s;
		if (i < 10)
			s = "D:/sandbox/fiberSimulation/hairs/scaledHairs/frame0000" + std::to_string(i) + "_scaled.txt";
		else 
			s = "D:/sandbox/fiberSimulation/hairs/scaledHairs/frame000" + std::to_string(i) + "_scaled.txt";
		const char* yarnfile1 = s.c_str();
		std::ifstream fin(yarnfile1);
		std::cout << yarnfile1 << std::endl;
		assert(fin.is_open() && "file wasn't found!\n");

		Fiber::Yarn yarn_tmp;
		const char* centerYarn1 = "centerYarn_ref.txt";
		yarn_tmp.yarnCenter(yarnfile1, centerYarn1);
		std::vector<yarnIntersect2D> pnts;
		CrossSection cs2(yarnfile1, centerYarn1, ply_num, n, 100, pnts);

		std::string t = "../data/allCrossSection2D_deformed_frame" + std::to_string(i) + ".txt";
		const char* deformed = t.c_str();
		plotIntersections(pnts, deformed);
	}

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
	//ply_centers_test();

	//extractNormals();

	//shapeMatch_test();
	//yarnShapeMatch_test();

	//writeNormals();

	render1fiber();
#endif

//	std::system("pause"); //add breakpoint instead


 	return 0;
}

void plotIntersections(const std::vector<yarnIntersect2D> &its, const char* filename) {
	FILE *fout;
	// write the plycenters
	if (fopen_s(&fout, filename, "wt") == 0) {
		const int ignorPlanes = 0.1 * its.size(); // crop the first and last 10% of the yarn
														
		fprintf_s(fout, "plane_num: %d \n", its.size() - 2 * ignorPlanes);
		fprintf_s(fout, "ply_num: %d \n", its[0].size());
		fprintf_s(fout, "\n");

		for (int i = ignorPlanes; i < its.size() - ignorPlanes; ++i) { //number of planes
			fprintf_s(fout, "index plane : %d \n", i - ignorPlanes);
			for (int p = 0; p < its[i].size(); ++p) { //number of plys
				fprintf_s(fout, "ply_fiber_num: %d \n", its[i][p].size());
				vec2f plyCenter(0.f);
				for (int j = 0; j < its[i][p].size(); ++j) { //number of intersections
					plyCenter += its[i][p][j];
				}
				plyCenter /= its[i][p].size();
				fprintf_s(fout, "plyCenter: %.4lf %.4lf \n", plyCenter.x, plyCenter.y);

				for (int j = 0; j < its[i][p].size(); ++j) { //number of intersections
					fprintf_s(fout, "%.4f %.4f \n", its[i][p][j].x, its[i][p][j].y);
				}
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
}

void plotIntersectionsDeformed(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Ellipse> &ellipses, const std::vector<float> &all_theta_R, const char* filename) {
	FILE *fout;
	// write the plycenters
	if (fopen_s(&fout, filename, "wt") == 0) {
		const int ignorPlanes = 0.1 * its.size(); // crop the first and last 10% of the yarn

		fprintf_s(fout, "plane_num: %d \n", its.size() - 2 * ignorPlanes);
		fprintf_s(fout, "ply_num: %d \n", its[0].size());
		fprintf_s(fout, "\n");

		Eigen::Matrix2f R, S, V, sigma, transf;
		for (int i = ignorPlanes; i < its.size() - ignorPlanes; ++i) { //number of planes
			sigma << ellipses[i].longR, 0, 0, ellipses[i].shortR;
			V << cos(ellipses[i].angle), -sin(ellipses[i].angle), sin(ellipses[i].angle), cos(ellipses[i].angle);
			S = V*sigma*V.transpose();
			R << cos(all_theta_R[i]), -sin(all_theta_R[i]), sin(all_theta_R[i]), cos(all_theta_R[i]);
			transf = R *S;

			fprintf_s(fout, "index plane : %d \n", i - ignorPlanes);
			for (int p = 0; p < its[i].size(); ++p) { //number of plys
				fprintf_s(fout, "ply_fiber_num: %d \n", its[i][p].size());
				vec2f plyCenter(0.f);
				for (int j = 0; j < its[i][p].size(); ++j) { //number of intersections
					plyCenter += its[i][p][j];
				}
				plyCenter /= its[i][p].size();
				Eigen::VectorXf ref(2, 1);
				ref << plyCenter.x, plyCenter.y;
				Eigen::VectorXf deformedRef = transf*ref;
				fprintf_s(fout, "plyCenter: %.6lf %.6lf \n", deformedRef[0], deformedRef[1]);

				for (int j = 0; j < its[i][p].size(); ++j) { //number of intersections
					Eigen::VectorXf ref(2, 1);
					ref << its[i][p][j].x, its[i][p][j].y;
					Eigen::VectorXf deformedRef = transf*ref;
					fprintf_s(fout, "%.6f %.6f \n", deformedRef[0], deformedRef[1]);
				}
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
}


void fittingCompress(CrossSection & cs, std::vector<yarnIntersect2D> &allPlaneIntersect, const char* compressFile) {

	//transfer from e1-e2 to x-y plane
	std::vector<yarnIntersect2D> xy_Its;
	cs.transferLocal2XY(allPlaneIntersect, xy_Its);

	//find the fitted ellipse parameters
	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(xy_Its, ellipses);

	//regularize the ellipse params
	std::vector<Ellipse> regulEllipses;
	//cs.regularizeEllipses(ellipses, regulEllipses, 20);
	//cs.regularizeEllipses(regulEllipses, ellipses, 20);
	//cs.regularizeEllipses(ellipses, regulEllipses, 20);
	//cs.regularizeEllipses(regulEllipses, ellipses, 20);

	//ellipses = regulEllipses;

	//update R1, R2, theta
	/*************** TO DO ***************/

	//do the dynamic programming optimization
	/*************** TO DO ***************/
	//optimizeEllipses()

	//parameterize the ellipse params
	//std::vector<Ellipse> param_ellipses;
	//cs.parameterizeEllipses(ellipses, param_ellipses);

	FILE *fout;
	if (fopen_s(&fout, compressFile, "wt") == 0) {
		fprintf_s(fout, "%d \n", ellipses.size());
		for (int i = 0; i < ellipses.size(); ++i) {
			fprintf_s(fout, "%.4f %.4f %.4f \n", ellipses[i].longR, ellipses[i].shortR, ellipses[i].angle);
		}
		fclose(fout);
	}

	FILE *fout1;
	//write ellipses to file for testing
	if (fopen_s(&fout1, "../data/orientation.txt", "wt") == 0) {
		const int ignorPlanes = 0.1 * ellipses.size(); // crop the first and last 10% of the yarn
		for (int i = ignorPlanes; i < ellipses.size() - ignorPlanes; ++i) {
			fprintf_s(fout1, "%.4f %.4f \n", 0.f, 0.f);
			fprintf_s(fout1, "%.4f %.4f %.4f \n", ellipses[i].longR, ellipses[i].shortR, ellipses[i].angle);
			fprintf_s(fout1, "\n");
		}
		fclose(fout1);
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