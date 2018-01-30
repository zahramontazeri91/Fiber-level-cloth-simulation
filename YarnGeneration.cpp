#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"
#include <string>

void phase1(const char* yarnfile1, const char* configfile, Fiber::Yarn &yarn, int skipFactor, int frame0, int frame1, int yarnNum, std::string &dataset) {
	std::cout << "*** Convert external force to local coordinate ***\n";

	for (int i = frame0; i < frame1; i++) {

		int f = i * skipFactor;

		HermiteCurve curve;
		int seg_subdiv = 100;

		for (int y = 0; y < yarnNum; ++y) {

			std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt"; //don't use upsampled centerline
			const char* curvefile = tmp7.c_str();
			std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";//don't use upsampled normals
			const char* normfile = tmp8.c_str();
			std::string tmp9 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(y) + "_world.txt";
			const char* physical_world = tmp9.c_str();
			std::string tmp10 = "input/" + dataset + "/physical_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local = tmp10.c_str();

			std::ifstream fin3(curvefile);
			assert(fin3.is_open() && "curvefile file wasn't found!\n");
			std::ifstream fin4(physical_world);
			assert(fin4.is_open() && "physical_world file wasn't found!\n");

			curve.init(curvefile, normfile, seg_subdiv);


			std::ifstream fin5(normfile);
			assert(fin5.is_open() && "normfile file wasn't found!\n");

			std::ifstream fin(physical_world);
			std::ofstream fout(physical_local);

			float S00, S01, S02, S10, S11, S12, S20, S21, S22;
			float A0, A1, A2;
			float B0, B1, B2;
			const int vrtx_num = yarn.getStepNum();
			for (int v = 0; v < vrtx_num; ++v) {
				fin >> S00 >> S01 >> S02 >> S10 >> S11 >> S12 >> S20 >> S21 >> S22
					>> A0 >> A1 >> A2
					>> B0 >> B1 >> B2;

				const double curveLength = curve.totalLength();
				float len = curveLength * (static_cast<double>(v) / static_cast<double>(vrtx_num - 1));
				const double t = curve.arcLengthInvApprox(len);

				Eigen::Vector3d ex, ey, ez;
				curve.getRotatedFrame(t, ex, ey, ez);

				/** local to world **/
				Eigen::Matrix3f local, world;
				world << S00, S01, S02,
					S10, S11, S12,
					S20, S21, S22;

				Eigen::Matrix3f M;
				M << ex[0], ex[1], ex[2],
					ey[0], ey[1], ey[2],
					ez[0], ez[1], ez[2];
				local = M*world*M.transpose();

				//write converted parameters
				fout << local(0, 0) << " " << local(0, 1) << " " << local(0, 2) << " " <<
					local(1, 0) << " " << local(1, 1) << " " << local(1, 2) << " " <<
					local(2, 0) << " " << local(2, 1) << " " << local(2, 2) << " ";

				Eigen::Vector3f localA, localB, worldA, worldB;
				worldA << A0, A1, A2;
				worldB << B0, B1, B2;
				localA = M*worldA;
				localB = M*worldB;
				fout << localA(0) << " " << localA(1) << " " << localA(2) << " " <<
					localB(0) << " " << localB(1) << " " << localB(2) << std::endl;
			}
			fout.close();
		}
	}

	std::cout << "*** Fitting phase ***\n";
	for (int i = frame0; i < frame1; i++) {
		int f = i * skipFactor;
		for (int y = 0; y < yarnNum; ++y) {

			//pipeline 2:
			//std::string tmp0 = "data/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + "_DEF.txt";
			//const char* yarnfile0 = tmp0.c_str();
			//std::string tmp2 = "input/" + dataset + "/deformGrad_" + std::to_string(cnt) + "_trans.txt";
			//const char* deformGrad = tmp2.c_str();

			std::string tmp1 = "data/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* yarnfile2 = tmp1.c_str();
			std::string tmp3 = "input/" + dataset + "/matrix_S_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S = tmp3.c_str();
			std::string tmp4 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			const char* curvefile = tmp4.c_str();
			std::string tmp5 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			const char* normfile = tmp5.c_str();

			std::ifstream fin1(yarnfile1);
			std::ifstream fin2(yarnfile2);
			//std::ifstream fin3(deformGrad);

			assert(fin1.is_open() && "reference-yarn file wasn't found!\n");
			assert(fin2.is_open() && "compressed-yarn file wasn't found!\n");
			//assert(fin3.is_open() && "deformGrad file wasn't found!\n");

			const int vrtx_num = yarn.getStepNum();

			//pipeline 2:
			//extractCompress_seg(configfile, yarnfile0, yarnfile0, deformGrad, compress_S,
			//curvefile, normfile, yarn.getPlyNum(), vrtx_num);
			extractCompress_seg(configfile, yarnfile1, yarnfile2, "noNeed.txt", compress_S,
				curvefile, normfile, yarn.getPlyNum(), vrtx_num);
			/*************************************************/
			std::string tmp6 = "output/" + dataset + "/genYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* outfile = tmp6.c_str();
			//// Procedural step
			yarn.simulate_ply();
			yarn.write_plys("test_ply.txt");
			const int K = yarn.getPlyNum();
			yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
			yarn.build("test_fly.txt", K);

			////pipeline 2:
			////yarn.compress_yarn3D(deformGrad, compress_S);

			yarn.compress_yarn_A(compress_S);
			yarn.curve_yarn(curvefile, normfile);
			yarn.write_yarn(outfile);
			///////*************************************************/
			std::string tmp7 = "output/" + dataset + "/genYarn_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* outfile_wo = tmp7.c_str();
			yarn.simulate_ply();
			yarn.write_plys("test_ply.txt");
			yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
			yarn.build("test_fly.txt", K);
			yarn.curve_yarn(curvefile, normfile);
			yarn.write_yarn(outfile_wo);
		}
	}
}
void phase2(const char* yarnfile1, const char* configfile, Fiber::Yarn &yarn, int skipFactor, int frame0, int frame1, int yarnNum, std::string &dataset) {
	std::cout << "*** Training phase ***\n";

	for (int i = frame0; i < frame1; i++) {

		int f = i * skipFactor;
		for (int y = 0; y < yarnNum; ++y) {

			std::string tmp6 = "input/" + dataset + "/NN/testY_NN_full_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			const char* curvefile = tmp7.c_str();
			std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			const char* normfile = tmp8.c_str();

			//std::string tmp9 = "input/" + dataset + "/deformGrad_" + std::to_string(f) + "_trans.txt";
			//const char* deformGrad = tmp9.c_str();

			std::cout << compress_S << std::endl;
			std::ifstream fin2(compress_S);
			assert(fin2.is_open() && "compress_S_NN file wasn't found!\n");
			std::ifstream fin3(curvefile);
			assert(fin3.is_open() && "curvefile file wasn't found!\n");
			std::ifstream fin4(normfile);
			assert(fin4.is_open() && "normfile file wasn't found!\n");

			std::string tmp3 = "output/" + dataset + "/genYarn_NN_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* outfile = tmp3.c_str();
			// Procedural step
			yarn.simulate_ply();
			yarn.write_plys("test_ply.txt");
			const int K = yarn.getPlyNum();
			yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
			yarn.build("test_fly.txt", K);
			//pipeline 2:
			//yarn.compress_yarn3D(deformGrad, compress_S);

			yarn.compress_yarn_A(compress_S);
			yarn.curve_yarn(curvefile, normfile);
			yarn.write_yarn(outfile);
			std::cout << outfile << std::endl;

			/*******  Validate NN by L2-norm ******/
			//std::string tmp4 = "output/" + dataset + "/genYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			//const char* yarnfile_proc = tmp4.c_str(); //proc yarn
			//std::ifstream fin6(yarnfile_proc);
			//assert(fin6.is_open() && "yarn_proc file wasn't found!\n");
			//Fiber::Yarn yarn_proc;
			//yarn_proc.parse(configfile);
			//yarn_proc.build(yarnfile_proc, yarn_proc.getPlyNum());

			//std::string tmp5 = "data/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			//const char* yarnfile_simul = tmp5.c_str();
			//std::ifstream fin5(yarnfile_simul);
			//assert(fin5.is_open() && "yarn_simul file wasn't found!\n");
			//Fiber::Yarn yarn_simul;
			//yarn_simul.parse(configfile);
			//yarn_simul.build(yarnfile_simul, yarn_simul.getPlyNum());

			//const int trimPercent = 0.15; // should match with building NN data
			//float L2;
			//yarn.L2norm_3D(yarn, yarn_proc, trimPercent, L2);
			//std::cout << "L2 is: " << L2 << std::endl;

		}
	}

}

int main(int argc, const char **argv) {

	const char* yarnfile1 = "genYarn_ref.txt";
	const char* configfile = "config_300.txt";
	std::ifstream fin0(configfile);
	assert(fin0.is_open() && "config file wasn't found!\n");
	Fiber::Yarn yarn;
	yarn.parse(configfile);
	yarn.yarn_simulate();
	yarn.write_yarn(yarnfile1);

	int yarnNum = 1;
	int skipFactor = 500;
	int frame0 = 8000 / skipFactor ;
	int frame1 = 15000 / skipFactor + 1;
	std::string dataset = "spacing1.0x_00011_straight";

	int phase = 8;

	switch (phase) {
	case 1: {
		std::cout << "*** Convert external force to local coordinate ***\n";

		for (int i = frame0; i < frame1; i++) {

			int f = i * skipFactor;

			HermiteCurve curve;
			int seg_subdiv = 100;

			for (int y = 0; y < yarnNum; ++y) {

				std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt"; //don't use upsampled centerline
				const char* curvefile = tmp7.c_str();
				std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";//don't use upsampled normals
				const char* normfile = tmp8.c_str();
				std::string tmp9 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(y) + "_world.txt";
				const char* physical_world = tmp9.c_str();
				std::string tmp10 = "input/" + dataset + "/physical_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				const char* physical_local = tmp10.c_str();

				std::ifstream fin3(curvefile);
				assert(fin3.is_open() && "curvefile file wasn't found!\n");
				std::ifstream fin4(physical_world);
				assert(fin4.is_open() && "physical_world file wasn't found!\n");

				curve.init(curvefile, normfile, seg_subdiv);


				std::ifstream fin5(normfile);
				assert(fin5.is_open() && "normfile file wasn't found!\n");

				std::ifstream fin(physical_world);
				std::ofstream fout(physical_local);

				float S00, S01, S02, S10, S11, S12, S20, S21, S22;
				float A0, A1, A2;
				float B0, B1, B2;
				const int vrtx_num = yarn.getStepNum();
				for (int v = 0; v < vrtx_num; ++v) {
					fin >> S00 >> S01 >> S02 >> S10 >> S11 >> S12 >> S20 >> S21 >> S22
						>> A0 >> A1 >> A2
						>> B0 >> B1 >> B2;

					const double curveLength = curve.totalLength();
					float len = curveLength * (static_cast<double>(v) / static_cast<double>(vrtx_num - 1));
					const double t = curve.arcLengthInvApprox(len);

					Eigen::Vector3d ex, ey, ez;
					curve.getRotatedFrame(t, ex, ey, ez);

					/** local to world **/
					Eigen::Matrix3f local, world;
					world << S00, S01, S02,
						S10, S11, S12,
						S20, S21, S22;

					Eigen::Matrix3f M;
					M << ex[0], ex[1], ex[2],
						ey[0], ey[1], ey[2],
						ez[0], ez[1], ez[2];
					local = M*world*M.transpose();

					//write converted parameters
					fout << local(0, 0) << " " << local(0, 1) << " " << local(0, 2) << " " <<
						local(1, 0) << " " << local(1, 1) << " " << local(1, 2) << " " <<
						local(2, 0) << " " << local(2, 1) << " " << local(2, 2) << " ";

					Eigen::Vector3f localA, localB, worldA, worldB;
					worldA << A0, A1, A2;
					worldB << B0, B1, B2;
					localA = M*worldA;
					localB = M*worldB;
					fout << localA(0) << " " << localA(1) << " " << localA(2) << " " <<
						localB(0) << " " << localB(1) << " " << localB(2) << std::endl;
				}
				fout.close();
			}
		}

		std::cout << "*** Fitting phase ***\n";
		for (int i = frame0; i < frame1; i++) {
			int f = i * skipFactor;
			for (int y = 0; y < yarnNum; ++y) {

				//pipeline 2:
				//std::string tmp0 = "data/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + "_DEF.txt";
				//const char* yarnfile0 = tmp0.c_str();
				//std::string tmp2 = "input/" + dataset + "/deformGrad_" + std::to_string(cnt) + "_trans.txt";
				//const char* deformGrad = tmp2.c_str();


				std::string tmp1 = "data/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				const char* yarnfile2 = tmp1.c_str();
				std::string tmp3 = "input/" + dataset + "/matrix_S_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				const char* compress_S = tmp3.c_str();
				std::string tmp4 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
				const char* curvefile = tmp4.c_str();
				std::string tmp5 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
				const char* normfile = tmp5.c_str();

				std::ifstream fin1(yarnfile1);
				std::ifstream fin2(yarnfile2);
				//std::ifstream fin3(deformGrad);

				assert(fin1.is_open() && "reference-yarn file wasn't found!\n");
				assert(fin2.is_open() && "compressed-yarn file wasn't found!\n");
				//assert(fin3.is_open() && "deformGrad file wasn't found!\n");

				const int vrtx_num = yarn.getStepNum();

				//pipeline 2:
				//extractCompress_seg(configfile, yarnfile0, yarnfile0, deformGrad, compress_S,
				//curvefile, normfile, yarn.getPlyNum(), vrtx_num);
				extractCompress_seg(configfile, yarnfile1, yarnfile2, "noNeed.txt", compress_S,
					curvefile, normfile, yarn.getPlyNum(), vrtx_num);
				/*************************************************/
				std::string tmp6 = "output/" + dataset + "/genYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				const char* outfile = tmp6.c_str();
				//// Procedural step
				yarn.simulate_ply();
				yarn.write_plys("test_ply.txt");
				const int K = yarn.getPlyNum();
				yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
				yarn.build("test_fly.txt", K);

				////pipeline 2:
				////yarn.compress_yarn3D(deformGrad, compress_S);

				yarn.compress_yarn_A(compress_S);
				yarn.curve_yarn(curvefile, normfile);
				yarn.write_yarn(outfile);
				///////*************************************************/
				std::string tmp7 = "output/" + dataset + "/genYarn_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				const char* outfile_wo = tmp7.c_str();
				yarn.simulate_ply();
				yarn.write_plys("test_ply.txt");
				yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
				yarn.build("test_fly.txt", K);
				yarn.curve_yarn(curvefile, normfile);
				yarn.write_yarn(outfile_wo);
			}
		}
		break;
	}
	case 2: {
		std::cout << "*** Training phase ***\n";

		for (int i = frame0; i < frame1; i++) {

			int f = i * skipFactor;
			for (int y = 0; y < yarnNum; ++y) {

				std::string tmp6 = "input/" + dataset + "/NN/testY_NN_full_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				const char* compress_S = tmp6.c_str();
				std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
				const char* curvefile = tmp7.c_str();
				std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
				const char* normfile = tmp8.c_str();

				//std::string tmp9 = "input/" + dataset + "/deformGrad_" + std::to_string(f) + "_trans.txt";
				//const char* deformGrad = tmp9.c_str();
				std::cout << compress_S << std::endl;
				std::ifstream fin2(compress_S);
				assert(fin2.is_open() && "compress_S_NN file wasn't found!\n");
				std::ifstream fin3(curvefile);
				assert(fin3.is_open() && "curvefile file wasn't found!\n");
				std::ifstream fin4(normfile);
				assert(fin4.is_open() && "normfile file wasn't found!\n");

				std::string tmp3 = "output/" + dataset + "/genYarn_NN_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				const char* outfile = tmp3.c_str();
				// Procedural step
				yarn.simulate_ply();
				yarn.write_plys("test_ply.txt");
				const int K = yarn.getPlyNum();
				yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
				yarn.build("test_fly.txt", K);
				//pipeline 2:
				//yarn.compress_yarn3D(deformGrad, compress_S);

				yarn.compress_yarn_A(compress_S);
				yarn.curve_yarn(curvefile, normfile);
				yarn.write_yarn(outfile);
				std::cout << outfile << std::endl;

				///*******  Validate NN by L2-norm ******/
				//std::string tmp4 = "output/" + dataset + "/genYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				//const char* yarnfile_proc = tmp4.c_str(); //proc yarn
				//std::ifstream fin6(yarnfile_proc);
				//assert(fin6.is_open() && "yarn_proc file wasn't found!\n");
				//Fiber::Yarn yarn_proc;
				//yarn_proc.parse(configfile);
				//yarn_proc.build(yarnfile_proc, yarn_proc.getPlyNum());

				//std::string tmp5 = "data/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				//const char* yarnfile_simul = tmp5.c_str(); 
				//std::ifstream fin5(yarnfile_simul);
				//assert(fin5.is_open() && "yarn_simul file wasn't found!\n");
				//Fiber::Yarn yarn_simul;
				//yarn_simul.parse(configfile);
				//yarn_simul.build(yarnfile_simul, yarn_simul.getPlyNum());

				//const int trimPercent = 0.15; // should match with building NN data
				//float L2;
				//yarn.L2norm_3D(yarn, yarn_proc, trimPercent, L2);
				//std::cout << "L2 is: " << L2 << std::endl;

			}
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
						fprintf_s(foutS, "%.6f %.6f %.6f ", (all_mat_S[i][j](0, 0) + all_mat_S[i + 1][j](0, 0)) / 2.f,
							(all_mat_S[i][j](1, 1) + all_mat_S[i + 1][j](1, 1)) / 2.f,
							(all_mat_S[i][j](0, 1) + all_mat_S[i + 1][j](0, 1)) / 2.f);
					}
					fprintf_s(foutS, "\n");
				}
				int i = all_mat_S.size() - 1;
				for (int j = 0; j < all_mat_S[i].size(); ++j) {
					fprintf_s(foutS, "%.6f %.6f %.6f", all_mat_S[i][j](0, 0) / 2.f,
						all_mat_S[i][j](1, 1) / 2.f,
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
						fprintf_s(foutT, "%.6f %.6f ", all_T[i][j](0, 0), all_T[i][j](1, 0));
					fprintf_s(foutT, "\n");
				}
				fclose(foutT);
			}



			/**************************************************/
			//std::string tmp6 = "output/" + dataset + "/genYarn_" + std::to_string(i * 5) + ".txt";
			//const char* outfile = tmp6.c_str();
			////yarn.compress_yarn(compress_R, compress_S);
			//yarn.compress_yarn(all_mat_S, all_theta_R, all_T);
			//yarn.curve_yarn(curvefile, normfile);
			//yarn.write_yarn(outfile);

			std::string tmp7 = "output/" + dataset + "/genYarn_wo_" + std::to_string(i * 5) + ".txt";
			const char* outfile_wo = tmp7.c_str();
			yarn.curve_yarn(curvefile, normfile);
			yarn.write_yarn(outfile_wo);
		}

		break;
	}
	case 4: {
		std::cout << "***Mapping yarn to a curve *** \n";
		dataset = "spacing0.5x_00011";
		std::string datatype = "fixed_teeth";
		yarnNum = 1;

		//Generate a yarn mapped to a given curve
		const char* configFile = "config_300.txt";
		std::ifstream fin1(configFile);
		assert(fin1.is_open() && "config file wasn't found!\n");

		Fiber::Yarn yarn0;
		yarn0.parse(configFile);
		yarn0.setStepNum(121);
		
		for (int y = 0; y < yarnNum; ++y) {
			std::string tmp1 = "../../dataSets/" + datatype + "/test/" + dataset + "/curves/curve_0_" + std::to_string(y) + ".txt";
			const char* curvefile = tmp1.c_str();
			//const char* normfile = "normYarn.txt";
			std::ifstream fin2(curvefile);
			assert(fin2.is_open() && "curve file wasn't found!\n");

			// Procedural step
			std::string ind = "";
			if (y < 10)
				ind = "0" + std::to_string(y);
			else 
				ind = std::to_string(y);
			std::string tmp2 = "../../dataSets/" + datatype + "/test/" + dataset + "/fiber/frame_0000000fiber_" + ind + ".obj";
			const char* simul_frame0 = tmp2.c_str();
			yarn0.yarn_simulate();
			yarn0.curve_yarn(curvefile);
			yarn0.write_yarn_obj(simul_frame0);
		}
		break;
	}
	case 5: {
		std::cout << "*** Convert external force to local coordinate ***\n";
		for (int i = frame0; i < frame1; i++) {
			int f = i * skipFactor;
			HermiteCurve curve;
			int seg_subdiv = 100;
			for (int y = 0; y < yarnNum; ++y) {

				std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt"; //don't use upsampled centerline
				const char* curvefile = tmp7.c_str();
				std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";//don't use upsampled normals
				const char* normfile = tmp8.c_str();
				std::string tmp9 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(y) + "_world.txt";
				const char* physical_world = tmp9.c_str();
				std::string tmp10 = "input/" + dataset + "/physical_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
				const char* physical_local = tmp10.c_str();

				std::ifstream fin3(curvefile);
				assert(fin3.is_open() && "curvefile file wasn't found!\n");
				std::ifstream fin4(physical_world);
				assert(fin4.is_open() && "physical_world file wasn't found!\n");

				curve.init(curvefile, normfile, seg_subdiv);


				std::ifstream fin5(normfile);
				assert(fin5.is_open() && "normfile file wasn't found!\n");

				std::ifstream fin(physical_world);
				std::ofstream fout(physical_local);

				float S00, S01, S02, S10, S11, S12, S20, S21, S22;
				float A0, A1, A2;
				float B0, B1, B2;
				const int vrtx_num = yarn.getStepNum();
				for (int v = 0; v < vrtx_num; ++v) {
					fin >> S00 >> S01 >> S02 >> S10 >> S11 >> S12 >> S20 >> S21 >> S22
						>> A0 >> A1 >> A2
						>> B0 >> B1 >> B2;

					const double curveLength = curve.totalLength();
					float len = curveLength * (static_cast<double>(v) / static_cast<double>(vrtx_num - 1));
					const double t = curve.arcLengthInvApprox(len);

					Eigen::Vector3d ex, ey, ez;
					curve.getRotatedFrame(t, ex, ey, ez);

					/** local to world **/
					Eigen::Matrix3f local, world;
					world << S00, S01, S02,
						S10, S11, S12,
						S20, S21, S22;

					Eigen::Matrix3f M;
					M << ex[0], ex[1], ex[2],
						ey[0], ey[1], ey[2],
						ez[0], ez[1], ez[2];
					local = M*world*M.transpose();

					//write converted parameters
					fout << local(0, 0) << " " << local(0, 1) << " " << local(0, 2) << " " <<
						local(1, 0) << " " << local(1, 1) << " " << local(1, 2) << " " <<
						local(2, 0) << " " << local(2, 1) << " " << local(2, 2) << " ";

					Eigen::Vector3f localA, localB, worldA, worldB;
					worldA << A0, A1, A2;
					worldB << B0, B1, B2;
					localA = M*worldA;
					localB = M*worldB;
					fout << localA(0) << " " << localA(1) << " " << localA(2) << " " <<
						localB(0) << " " << localB(1) << " " << localB(2) << std::endl;
				}
				fout.close();
			}
		}
		break;
	}
	case 6: {
		std::cout << " *** test flyaways *** ";

		yarn.simulate_ply();
		yarn.write_plys("test_ply.txt");
		const int K = yarn.getPlyNum();
		yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
		yarn.build("test_fly.txt", K);

		const char* compress_R = "matrix_R_175.txt";
		const char* compress_S = "matrix_S_175.txt";
		const char* curvefile = "centerYarn_175.txt";
		const char* normfile = "normYarn_175.txt";

		yarn.compress_yarn(compress_R, compress_S);
		yarn.curve_yarn(curvefile, normfile);
		yarn.write_yarn("test_fly.txt");
		break;
	}
	case 7: {
		std::cout << " *** new pipeline: generate yarns using yarn-level *** " << std::endl;
		int cnt = frame0 * skipFactor;
		for (int i = frame0; i < frame1; i++) {
			int f = i * skipFactor;
			for (int y = 0; y < yarnNum; ++y) {

				std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(cnt) + "_ds.txt";
				const char* curvefile = tmp1.c_str();
				std::string tmp2 = "input/" + dataset + "/normYarn_" + std::to_string(cnt) + "_ds.txt";
				const char* normfile = tmp2.c_str();
				std::string tmp3 = "input/" + dataset + "/deformGrad_" + std::to_string(cnt) + ".txt";
				const char* deformGrad_local = tmp3.c_str();

				std::ifstream fin1(curvefile);
				std::ifstream fin2(normfile);
				std::ifstream fin3(deformGrad_local);

				assert(fin1.is_open() && "curvefile file wasn't found!\n");
				assert(fin2.is_open() && "compressed-yarn file wasn't found!\n");
				assert(fin3.is_open() && "physical_local file wasn't found!\n"); 


				std::string tmp4 = "data/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + "_DEF.txt";
				const char* outfile_wo = tmp4.c_str();

				yarn.simulate_ply();
				yarn.write_plys("test_ply.txt");
				const int K = yarn.getPlyNum();
				yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
				yarn.build("test_fly.txt", K);

				//std::string tmp5 = "data/" + dataset + "/simul_frame_" + std::to_string(0) + "_" + std::to_string(y) + "_4x.txt";
				//const char* yarnfile = tmp5.c_str();
				//yarn.parse(configfile);
				//yarn.build(yarnfile, 2);

				
				yarn.compress_yarn3D(deformGrad_local);
				yarn.write_yarn(outfile_wo);

				cnt += skipFactor;
			}
		}
		break;
	}
	case 8: {
		/**************** RUN ALL ****************/
		frame0 = 8000 / skipFactor ;
		
		dataset = "spacing0.5x";
		frame1 = 0 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing1.0x";
		frame1 = 0 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing1.5x";
		frame1 = 15000 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);

		dataset = "spacing0.5x_00011";
		frame1 = 16000 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing0.5x_10100";
		frame1 = 15000 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing0.5x_11110";
		frame1 = 15000 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing1.0x_00011";
		frame1 = 17000 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing1.0x_10100";
		frame1 = 15500 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing1.0x_11110";
		frame1 = 16000 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing1.5x_00011";
		frame1 = 17500 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing1.5x_10100";
		frame1 = 16000 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		dataset = "spacing1.5x_11110";
		frame1 = 16500 / skipFactor + 1;
		phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);

		/********************************/
		break;
	}
	case 9: {
		/**************** RUN ALL ****************/
		//frame0 = 5 / skipFactor + 1;
		//dataset = "spacing1.0x_p1";
		//frame1 = 30 / skipFactor + 1;
		//phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);

		frame0 = 140 / skipFactor + 1;
		dataset = "spacing0.5x_00011";
		frame1 = 160 / skipFactor + 1;
		phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		//dataset = "spacing1.5x_10100";
		//frame1 = 160 / skipFactor + 1;
		//phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		//dataset = "spacing1.5x_11110";
		//frame1 = 165 / skipFactor + 1;
		//phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		//dataset = "spacing1.0x_00011";
		//frame1 = 170 / skipFactor + 1;
		//phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		//dataset = "spacing1.0x_10100";
		//frame1 = 155 / skipFactor + 1;
		//phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		//dataset = "spacing1.0x_11110";
		//frame1 = 160 / skipFactor + 1;
		//phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
		/********************************/
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