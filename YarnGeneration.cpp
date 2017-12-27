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
	
	int frame0 = 31;
	int frame1 = 32;
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
				
				//return 0;

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
			std::cout << "*** Per ply shape matching  ***\n";
			for (int i = frame0; i < frame1; i++) {

				// Procedural step
				yarn.yarn_simulate();

				std::string tmp1 = "data/" + dataset + "/simul_frame_" + std::to_string(i * 5) + ".txt";
				const char* yarnfile2 = tmp1.c_str();
				std::string tmp2 = "input/" + dataset + "/matrix_R_" + std::to_string(i * 5) + ".txt";
				const char* compress_R = tmp2.c_str();
				std::string tmp3 = "input/" + dataset + "/matrix_S_" + std::to_string(i * 5) + ".txt";
				const char* compress_S = tmp3.c_str();

				std::string tmp = "input/" + dataset + "/matrix_T_" + std::to_string(i * 5) + ".txt";
				const char* compress_T = tmp.c_str();

				std::string tmp4 = "input/" + dataset + "/centerYarn_" + std::to_string(i * 5) + ".txt";
				const char* curvefile = tmp4.c_str();
				std::string tmp5 = "input/" + dataset + "/normYarn_" + std::to_string(i * 5) + ".txt";
				const char* normfile = tmp5.c_str();

				std::ifstream fin1(yarnfile1);
				std::ifstream fin2(yarnfile2);
				assert(fin1.is_open() && "reference-yarn file wasn't found!\n");
				assert(fin2.is_open() && "compressed-yarn file wasn't found!\n");

				/*************************************************/
				std::vector<yarnIntersect2D> pnts_ref;
				yarn.yarn2crossSections(pnts_ref);

				Fiber::Yarn yarn_tmp;
				yarn_tmp.yarnCenter(yarnfile2, curvefile);
				std::vector<yarnIntersect2D> pnts_trans;
				CrossSection cs2(yarnfile2, curvefile, normfile, yarn.getPlyNum(), yarn.getStepNum(), 100, pnts_trans, false);

				std::vector<std::vector<Eigen::MatrixXf>> all_mat_S;
				std::vector<std::vector<float>> all_theta_R;
				std::vector<std::vector<Eigen::MatrixXf>> all_T;
				cs2.yarnShapeMatches(pnts_trans, pnts_ref, all_mat_S, all_theta_R, all_T);

				FILE *foutR;
				if (fopen_s(&foutR, compress_R, "wt") == 0) {
					fprintf_s(foutR, "%d \n", all_theta_R.size()); //number of planes
					for (int i = 0; i < all_theta_R.size(); ++i) {
						//fprintf_s(foutR, "%d \n", all_theta_R[i].size()); //number of plys
						for (int j = 0; j < all_theta_R[i].size(); ++j) 
							fprintf_s(foutR, "%.6f ", all_theta_R[i][j]);
						fprintf_s(foutR, "\n");
					}
					fclose(foutR);
				}

				FILE *foutS;
				// write S-matrix for each segment not vertex 
				if (fopen_s(&foutS, compress_S, "wt") == 0) {
					//fprintf_s(foutS, "%d \n", all_mat_S.size());
					for (int i = 0; i < all_mat_S.size() - 1; ++i) { //number of planes
						for (int j = 0; j < all_mat_S[i].size(); ++j) { //number of plys
							fprintf_s(foutS, "%.6f %.6f %.6f ", (all_mat_S[i][j](0,0) + all_mat_S[i + 1][j](0, 0)) / 2.f,
								(all_mat_S[i][j](1,1) + all_mat_S[i + 1][j](1,1)) / 2.f,
								(all_mat_S[i][j](0, 1) + all_mat_S[i + 1][j](0, 1)) / 2.f);
						}
						fprintf_s(foutS, "\n");
					}
					int i = all_mat_S.size() - 1;
					for (int j = 0; j < all_mat_S[i].size(); ++j) {
						fprintf_s(foutS, "%.6f %.6f %.6f", all_mat_S[i][j](0, 0) / 2.f,
							all_mat_S[i][j](1,1) / 2.f,
							all_mat_S[i][j](0, 1) / 2.f);
					}
					fprintf_s(foutS, "\n");
					fclose(foutS);
				}

				FILE *foutT;
				if (fopen_s(&foutT, compress_T, "wt") == 0) {
					fprintf_s(foutT, "%d \n", all_T.size()); //number of planes
					for (int i = 0; i < all_T.size(); ++i) {
						//fprintf_s(foutT, "%d \n", all_theta_R[i].size()); //number of plys
						for (int j = 0; j < all_T[i].size(); ++j)
							fprintf_s(foutT, "%.6f %.6f ", all_T[i][j](0,0), all_T[i][j](1, 0));
						fprintf_s(foutT, "\n");
					}
					fclose(foutT);
				}



				/**************************************************/
				std::string tmp6 = "output/" + dataset + "/genYarn_" + std::to_string(i * 5) + ".txt";
				const char* outfile = tmp6.c_str();

				//yarn.compress_yarn(compress_R, compress_S);
				yarn.compress_yarn(all_mat_S, all_theta_R, all_T);
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
		case 4: {
			std::cout << "*** Generation phase (Chang use) *** \n";
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