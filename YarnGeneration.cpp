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
		int seg_subdiv = 10;
		for (int y = 0; y < yarnNum; ++y) {

			//std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt"; //don't use upsampled centerline
			//const char* curvefile = tmp7.c_str();
			//std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";//don't use upsampled normals
			//const char* normfile = tmp8.c_str();
			std::string tmp9 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(y) + "_world.txt";
			const char* physical_world = tmp9.c_str();
			std::string tmp10 = "input/" + dataset + "/physical_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local = tmp10.c_str();

			std::string tmp1 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* normfile_us = tmp1.c_str();
			std::string tmp2 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* curvefile_us = tmp2.c_str();

			std::string tmp5 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* twistfile = tmp5.c_str();

			std::ifstream fin2(curvefile_us);
			assert(fin2.is_open() && "curvefile_us file wasn't found!\n");
			//std::ifstream fin3(curvefile);
			//assert(fin3.is_open() && "curvefile file wasn't found!\n");
			std::ifstream fin4(physical_world);
			assert(fin4.is_open() && "physical_world file wasn't found!\n");

			curve.init_norm(curvefile_us, normfile_us, seg_subdiv); /******* Change upsample to downsample *******/
			std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;
			curve.assign(all_pts, all_tang, all_norm);
			//curve.assign_twist(twistfile, all_pts, all_tang, all_norm, 2);
			assert(all_pts.size() == yarn.getStepNum());
			//curve.init_principleNormal(curvefile, normfile, seg_subdiv);


			//std::ifstream fin5(normfile);
			//assert(fin5.is_open() && "normfile file wasn't found!\n");  

			std::ifstream fin(physical_world);
			std::ofstream fout(physical_local);

			float S00, S01, S02, S10, S11, S12, S20, S21, S22;
			float A0, A1, A2;
			float B0, B1, B2;
			const int vrtx_num = yarn.getStepNum();

			for (int v = 0; v < vrtx_num; ++v) {
				fin >> S00 >> S01 >> S02 >> S10 >> S11 >> S12 >> S20 >> S21 >> S22
					//>> A0 >> A1 >> A2
					//>> B0 >> B1 >> B2
					;

				/*
				const double curveLength = curve.totalLength();
				float len = curveLength * (static_cast<double>(v) / static_cast<double>(vrtx_num - 1));
				const double t = curve.arcLengthInvApprox(len);
				Eigen::Vector3d ex, ey, ez;
				curve.getRotatedFrame(t, ex, ey, ez); */


				Eigen::Vector3d ez = all_tang[v];
				Eigen::Vector3d ey = all_norm[v];
				Eigen::Vector3d ex = ez.cross(ey);

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

				/// uncomment for when having the internal forces
				//Eigen::Vector3f localA, localB, worldA, worldB;
				//worldA << A0, A1, A2;
				//worldB << B0, B1, B2;
				//localA = M*worldA;
				//localB = M*worldB;
				//fout << localA(0) << " " << localA(1) << " " << localA(2) << " " <<
				//localB(0) << " " << localB(1) << " " << localB(2) ;

				fout << std::endl;
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
			//std::string tmp4 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			//const char* curvefile = tmp4.c_str();
			//std::string tmp5 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			//const char* normfile = tmp5.c_str();

			std::string tmp6 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* normfile_us = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* curvefile_us = tmp7.c_str();

			std::string tmp5 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* twistfile = tmp5.c_str();

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
			const int upsample = 2;
			extractCompress_seg(configfile, yarnfile1, yarnfile2, "noNeed.txt", compress_S,
				curvefile_us, normfile_us, twistfile, 2, yarn.getPlyNum(), vrtx_num);
			/*************************************************/
			std::string tmp8 = "output/" + dataset + "/genYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* outfile = tmp8.c_str();
			//// Procedural step
			yarn.simulate_ply();
			yarn.write_plys("test_ply.txt");
			const int K = yarn.getPlyNum();
			yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
			yarn.build("test_fly.txt", K);

			////pipeline 2:
			////yarn.compress_yarn3D(deformGrad, compress_S);

			yarn.compress_yarn_A(compress_S);
			yarn.curve_yarn(curvefile_us, normfile_us);
			yarn.write_yarn(outfile);

			///////*************************************************/
			//std::string tmp7 = "output/" + dataset + "/genYarn_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			//const char* outfile_wo = tmp7.c_str();
			//yarn.simulate_ply();
			//yarn.write_plys("test_ply.txt");
			//yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
			//yarn.build("test_fly.txt", K);
			//yarn.curve_yarn(curvefile, normfile);
			//yarn.write_yarn(outfile_wo);
		}
	}
}
void phase2(const char* yarnfile1, const char* configfile, Fiber::Yarn &yarn, int skipFactor, int frame0, int frame1, int yarnNum, std::string &dataset) {
	std::cout << "*** Training phase ***\n";

	for (int i = frame0; i < frame1; i++) {

		int f = i * skipFactor;
		for (int y = 0; y < yarnNum; ++y) {

			std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* curvefile_us = tmp1.c_str();
			std::string tmp2 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* normfile_us = tmp2.c_str();

			std::string tmp6 = "input/" + dataset + "/NN/testY_NN_full_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S = tmp6.c_str();
			//std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			//const char* curvefile = tmp7.c_str();
			//std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			//const char* normfile = tmp8.c_str();

			//std::string tmp9 = "input/" + dataset + "/deformGrad_" + std::to_string(f) + "_trans.txt";
			//const char* deformGrad = tmp9.c_str();
			std::cout << compress_S << std::endl;
			std::ifstream fin2(compress_S);
			assert(fin2.is_open() && "compress_S_NN file wasn't found!\n");
			std::ifstream fin3(curvefile_us);
			assert(fin3.is_open() && "curvefile file wasn't found!\n");
			std::ifstream fin4(normfile_us);
			assert(fin4.is_open() && "normfile file wasn't found!\n");

			///*******  write the yarn ******/
			std::string tmp3 = "output/" + dataset + "/genYarn_NN_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* outfile = tmp3.c_str();
			// Procedural step
			yarn.yarn_simulate();
			//pipeline 2:
			//yarn.compress_yarn3D(deformGrad, compress_S);
			yarn.compress_yarn_A(compress_S);
			yarn.curve_yarn(curvefile_us, normfile_us);
			yarn.write_yarn(outfile);
			//std::cout << outfile << std::endl;

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
void buildTraining(const char* curvefile_ds, const char* normfile_ds, const char* physical_world, const char* compress_S, Fiber::Yarn &yarn, const float trimPercent,
	const int window_size, const char* curvefile_us, const char* angles, const char* physical_local_window, const char* physical_local_window_test, const char* compress_S_window,
	std::ofstream &fout_trainX_all, std::ofstream &fout_trainY_all, int isTrain ) {

	std::ifstream fin1(curvefile_ds);
	assert(fin1.is_open() && "curvefile_ds file wasn't found!\n");
	std::ifstream fin3(physical_world);
	assert(fin3.is_open() && "physical_world file wasn't found!\n");
	std::ifstream fin4(compress_S);
	if (isTrain) assert(fin4.is_open() && "compress_S file wasn't found!\n");

	std::ifstream fin_dg(physical_world);
	std::ifstream fin_S(compress_S);
	std::ofstream fout_dg(physical_local_window);
	std::ofstream fout_S(compress_S_window);
	std::ofstream fout_cntr(curvefile_us);
	std::ofstream fout_angle(angles);

	std::ofstream fout_dg_test(physical_local_window_test);
	//debug: 
	std::ofstream fout_S_test("input/spacing1.0x_00011/NN/testY_17000_0.txt");
	const char* normfile_us = "input/spacing1.0x_00011/normYarn_17000_0_us.txt";
	std::ofstream fout_TNB_yarn("../data/TNB_yarn.txt");

	int seg_subdiv = 10;
	HermiteCurve fullCurve_ds;
	fullCurve_ds.init(curvefile_ds, normfile_ds, seg_subdiv);
	const int vrtx_num = yarn.getStepNum();
	// write up-sampled centerline so later we can crop a segment out of it
	const double fullCurveLength = fullCurve_ds.totalLength();
	std::vector<Eigen::Matrix3f> all_dg;
	std::vector<Eigen::Matrix2f> all_S;
	std::vector<Eigen::Vector3d> all_pnt;
	std::vector<Eigen::Vector3d> all_n, all_t;
	fout_cntr << vrtx_num << "\n";

#if 0
	//generate the up-sampled yarn
	for (int v = 0; v < vrtx_num; ++v) {
		float fullLen = fullCurveLength * (static_cast<double>(v) / static_cast<double>(vrtx_num - 1));
		const double t_fll = fullCurve_ds.arcLengthInvApprox(fullLen);
		Eigen::Vector3d pnt = fullCurve_ds.eval(t_fll);
		fout_cntr << pnt[0] << " " << pnt[1] << " " << pnt[2] << "\n";
	}
	fout_cntr.close();
	HermiteCurve fullCurve_us;
	fullCurve_us.init(curvefile_us, normfile_us, seg_subdiv);
#endif 

	for (int v = 0; v < vrtx_num; ++v) {
		float fullLen = fullCurveLength * (static_cast<double>(v) / static_cast<double>(vrtx_num - 1));
		const double t_fll = fullCurve_ds.arcLengthInvApprox(fullLen);

		//****store normals for all points
		Eigen::Vector3d p = fullCurve_ds.eval(t_fll);
		all_pnt.push_back(p);
		Eigen::Vector3d n = fullCurve_ds.evalNormal(t_fll);
		all_n.push_back(n);
		Eigen::Vector3d tg = fullCurve_ds.evalTangent(t_fll);
		all_t.push_back(tg);
		/* Note that normals don't exactly match with up-sampled curve because adding new vertices to curve changes its curvature a bit */

		//*****store dg for all points
		float dg00, dg01, dg02, dg10, dg11, dg12, dg20, dg21, dg22;
		fin_dg >> dg00 >> dg01 >> dg02 >> dg10 >> dg11 >> dg12 >> dg20 >> dg21 >> dg22;
		Eigen::Matrix3f world_dg;
		world_dg << dg00, dg01, dg02,
			dg10, dg11, dg12,
			dg20, dg21, dg22;
		all_dg.push_back(world_dg);

		//*******store shape-matching matrix
		if (isTrain) {
			//store S-matrix for all points
			float S00, S01, S10, S11;
			fin_S >> S00 >> S01 >> S10 >> S11;
			Eigen::Matrix2f S;
			S << S00, S01, S10, S11;
			all_S.push_back(S);
		}
	}

	const int ignorPlanes = trimPercent * vrtx_num; // crop the first and last #% of the yarn
	for (int w = ignorPlanes; w < (vrtx_num - window_size + 1) - ignorPlanes; w++) {
		//std::cout << w << std::endl;

		//define a curve segment 
		const int start = w;
		const int end = w + (window_size - 1);
		HermiteCurve curve;
		curve.init_seg(curvefile_ds, start, end, seg_subdiv);  ///NOTE HERE ###########
		const double curveLength = curve.totalLength();

		for (int v = 0; v < window_size; ++v) {
			const double curveLength = curve.totalLength();
			float len = curveLength * (static_cast<double>(v) / static_cast<double>(window_size - 1));
			const double t = curve.arcLengthInvApprox(len);


			Eigen::Vector3d ez = curve.evalTangent(t);
			Eigen::Vector3d ex = curve.evalNormal(t);
			Eigen::Vector3d ey = ez.cross(ex);

			/** local to world **/
			Eigen::Matrix3f local_dg, M;
			M << ex[0], ex[1], ex[2],
				ey[0], ey[1], ey[2],
				ez[0], ez[1], ez[2];
			const int indx = w + v;
			local_dg = M*all_dg[indx] * M.transpose();

			/* Debug: */
			//if (w == 1)
			//	std::cout << "\n --------------------------- \n" << curve.eval(t) << std::endl;
			//if (w == 101)
			//	std::cout << "\n +++++++++++++++++++++++++++ \n" << curve.eval(t) << std::endl;
			//if (w == 1)
			//	fout_dg << " \n --------------------------- \n ";
			//if (w == 51)
			//	fout_dg << " \n +++++++++++++++++++++++++++ \n ";

			//write converted parameters
			fout_dg << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
				local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
				local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
			//write converted parameters and accmulate for all frames 
			fout_trainX_all << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
				local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
				local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";

		}
		fout_dg << std::endl;
		fout_trainX_all << std::endl;

		const int v_full = ceil((start + end) / 2.0); //index for the full curve			
		Eigen::Vector3d n_full = all_n[v_full];

		const int v = ceil((end - start) / 2.0);
		float len = curveLength * (static_cast<double>(v) / static_cast<double>(window_size - 1));
		const double t = curve.arcLengthInvApprox(len);
		Eigen::Vector3d n = curve.evalNormal(t);
		assert(v_full == v + start && "index for full yarn must be equal to index segment added with starting vertex");

		// rotate the shape-matching matrix to align the new normal
		float angle = acos(n.dot(n_full));

#if 1
		Eigen::Vector3d cross = (n.cross(n_full)).normalized();
		angle = signbit(cross[2]) == signbit(curve.evalTangent(t)[2]) ? angle : -1.0*angle; //cross should be in same direction with tangent. If not, negate the angle
		//std::cout << signbit(cross[0]) << " " << cross[0] << std::endl << signbit(curve.evalTangent(t)[0]) << std::endl << curve.evalTangent(t)[0] << std::endl << std::endl;
#endif

		std::cout << cross << " \n " << curve.evalTangent(t) << std::endl << std::endl;

		Eigen::Matrix2f S_rot;
		if (isTrain) {
			Eigen::Matrix2f R, S;
			R << cos(angle), -sin(angle),
				sin(angle), cos(angle);
			S_rot = R*all_S[v_full] * R.transpose();
			//S_rot = all_S[v_full];
			fout_S << S_rot(0, 0) << " " << S_rot(0, 1) << " " << S_rot(1, 0) << " " << S_rot(1, 1) << "\n";
			//std::cout << S << std::endl << S_rot << std::endl << std::endl;
			fout_trainY_all << S_rot(0, 0) << " " << S_rot(0, 1) << " " << S_rot(1, 0) << " " << S_rot(1, 1) << "\n";
		}
		//fout_angle << angle << std::endl;

# if 0
		///****************************************************************/
		///***** augment the training data by rotating normals by 180 *****/
		for (int v = 0; v < window_size; ++v) {

			const double curveLength = curve.totalLength();
			float len = curveLength * (static_cast<double>(v) / static_cast<double>(window_size - 1));
			const double t = curve.arcLengthInvApprox(len);

			Eigen::Vector3d ez = curve.evalTangent(t);
			Eigen::Vector3d ex = -1.0 * curve.evalNormal(t);
			Eigen::Vector3d ey = ez.cross(ex);

			/** local to world **/
			Eigen::Matrix3f local_dg, M;
			M << ex[0], ex[1], ex[2],
				ey[0], ey[1], ey[2],
				ez[0], ez[1], ez[2];
			const int indx = w + v;
			local_dg = M*all_dg[indx] * M.transpose();

			//write converted parameters
			fout_dg << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
				local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
				local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
			//write converted parameters and accmulate for all frames 
			fout_trainX_all << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
				local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
				local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
		}
		fout_dg << std::endl;
		fout_trainX_all << std::endl;

		Eigen::Matrix2f S_rot_pi;
		if (isTrain) {
			Eigen::Matrix2f R, S;
			R << cos(angle), -sin(angle),
				sin(angle), cos(angle);
			S_rot_pi = R*all_S[v_full] * R.transpose();
			//S_rot_pi = all_S[v_full];
			fout_S << S_rot_pi(0, 0) << " " << S_rot_pi(0, 1) << " " << S_rot_pi(1, 0) << " " << S_rot_pi(1, 1) << "\n";
			fout_trainY_all << S_rot_pi(0, 0) << " " << S_rot_pi(0, 1) << " " << S_rot_pi(1, 0) << " " << S_rot_pi(1, 1) << "\n";

		}
		//fout_angle << angle + pi << std::endl;
#endif

		///****************************************************************/
		/*** write test data X ******/

		for (int v = 0; v < window_size; ++v) {

			const double curveLength = curve.totalLength();
			float len = curveLength * (static_cast<double>(v) / static_cast<double>(window_size - 1));
			const double t = curve.arcLengthInvApprox(len);

			Eigen::Vector3d ez, ey, ex;

			//if (angle < pi / 2.0) {
			//	ez = curve.evalTangent(t);
			//	ex = curve.evalNormal(t);
			//	ey = ez.cross(ex);
			//}
			//else {
			//	ez = curve.evalTangent(t);
			//	ex = -1.0 * curve.evalNormal(t);
			//	ey = ez.cross(ex);
			//}

			ez = curve.evalTangent(t);
			ex = curve.evalNormal(t);
			ey = ez.cross(ex);

			/** local to world **/
			Eigen::Matrix3f local_dg, M;
			M << ex[0], ex[1], ex[2],
				ey[0], ey[1], ey[2],
				ez[0], ez[1], ez[2];
			const int indx = w + v;
			local_dg = M*all_dg[indx] * M.transpose();

			//////////////////////////////////////
			//if (w == 1)
			//	fout_dg_test << " \n --------------------------- \n ";
			//if (w == 51)
			//	fout_dg_test << " \n +++++++++++++++++++++++++++ \n ";
			//write converted parameters
			fout_dg_test << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
				local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
				local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
		}
		//write testY
		//Eigen::Matrix2f tmp_S = angle < pi / 2.0 ? S_rot : S_rot_pi;
		//fout_S_test << tmp_S(0, 0) << " " << tmp_S(0, 1) << " " << tmp_S(1, 0) << " " << tmp_S(1, 1) << "\n";
		//float tmp = angle < pi/2.0 ? angle : angle + pi;
		//if (tmp != tmp) //dot product (must be 1) but might be larger than 1 and so acos return nan 
		//	tmp = n_full.dot(n) > 0 ? 0.0 : pi;
		//fout_angle << tmp << std::endl;

		/*****tmp ****/
		fout_S_test << S_rot(0, 0) << " " << S_rot(0, 1) << " " << S_rot(1, 0) << " " << S_rot(1, 1) << "\n";
		fout_angle << angle << std::endl;
		/****/

		fout_dg_test << std::endl;


	}
	
	fout_dg.close();
	fout_S.close();
	fout_angle.close();
	fout_S_test.close();
	fout_dg_test.close();
}
void buildTraning_all(Fiber::Yarn &yarn, int skipFactor, int frame0, int frame1, int yarnNum, std::string &dataset, const int window_size, const float trimPercent, const int isTrain ) {

	std::string tmp7 = "input/" + dataset + "/NN/trainX_all.txt";
	const char* all_trainX = tmp7.c_str();
	std::string tmp8 = "input/" + dataset + "/NN/trainY_all.txt";
	const char* all_trainY = tmp8.c_str();
	std::ofstream fout_trainX_all(all_trainX);
	std::ofstream fout_trainY_all(all_trainY);

	for (int i = frame0; i < frame1; i++) {

		int f = i * skipFactor;
		for (int y = 0; y < yarnNum; ++y) {

			std::string tmp0 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			const char* curvefile_ds = tmp0.c_str();
			std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* curvefile = tmp1.c_str();
			std::string tmp2 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
			const char* normfile_ds = tmp2.c_str();
			std::string tmp3 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(y) + "_world.txt";
			const char* physical_world = tmp3.c_str();
			std::string tmp4 = "input/" + dataset + "/NN/trainX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_window = tmp4.c_str();
			std::string tmp5 = "input/" + dataset + "/matrix_S_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S = tmp5.c_str();
			std::string tmp6 = "input/" + dataset + "/NN/trainY_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S_window = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/NN/angles_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* angles = tmp7.c_str();

			std::string tmp8 = "input/" + dataset + "/NN/testX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_window_test = tmp8.c_str();

			buildTraining(curvefile_ds, normfile_ds, physical_world, compress_S, yarn, trimPercent, window_size, curvefile, angles,
				physical_local_window, physical_local_window_test, compress_S_window, fout_trainX_all, fout_trainY_all, isTrain);
		}
	}

	fout_trainX_all.close();
	fout_trainY_all.close();
}

int main(int argc, const char **argv) {

	const char* yarnfile1 = "genYarn_ref.txt";
	const char* configfile = "config_300.txt";
	std::ifstream fin0(configfile);
	assert(fin0.is_open() && "config file wasn't found!\n");
	Fiber::Yarn yarn;
	yarn.parse(configfile);

	yarn.setStepNum(300);
	
	yarn.yarn_simulate();
	yarn.write_yarn(yarnfile1);

	int yarnNum = 1;
	int skipFactor = 1;
	int frame0 = 200 / skipFactor;
	int frame1 = 200 / skipFactor + 1;
	std::string dataset = "twist_only";
	//std::string dataset = "spacing1.0x_00011";

	int phase = 2;

	switch (phase) {
		case 0: {
			std::cout << "*** up-sample the curve ***\n";
			for (int i = frame0; i < frame1; i++) {
				int f = i * skipFactor;
				HermiteCurve curve_ds, curve_us;
				int seg_subdiv = 10;
				for (int y = 0; y < yarnNum; ++y) {

					std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt"; 
					const char* curvefile_ds = tmp1.c_str();
					std::string tmp2 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt"; 
					const char* curvefile_us = tmp2.c_str();
					std::string tmp3 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
					const char* normfile_ds = tmp3.c_str();
					std::string tmp4 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* normfile_us = tmp4.c_str();
					std::string tmp5 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* twistfile = tmp5.c_str();

					std::ofstream fout_cntr(curvefile_us);
					std::ofstream fout_norm(normfile_us);

					const int vrtx_num = yarn.getStepNum();
					int seg_subdiv = 10;
					HermiteCurve curve_ds;
					curve_ds.init(curvefile_ds, normfile_ds, seg_subdiv);
					std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;

					//curve_ds.assign(all_pts, all_tang, all_norm);
					curve_ds.assign_twist(twistfile, all_pts, all_tang, all_norm, 2);
					//curve_ds.assign_upsample(all_pts, all_tang, all_norm); //upsample by 2					
					//curve_ds.twistNormals(twistfile, all_tang, all_norm, all_norm_rot);

					// write up-sampled centerline so later we can crop a segment out of it
					fout_cntr << vrtx_num << std::endl;
					fout_norm << vrtx_num << std::endl;
					for (int v = 0; v < vrtx_num; ++v) {

						fout_cntr << all_pts[v][0] << " " << all_pts[v][1] << " " << all_pts[v][2] << "\n";
						fout_norm << all_norm[v][0] << " " << all_norm[v][1] << " " << all_norm[v][2] << "\n";
					}
					fout_cntr.close();
					fout_norm.close();

					std::ofstream fout_TNB("../data/TNB.txt");
					for (int i = 0; i < vrtx_num; i++) {
						fout_TNB << all_pts[i][0] << " " << all_pts[i][1] << " " << all_pts[i][2] << " " <<
							all_tang[i][0] << " " << all_tang[i][1] << " " << all_tang[i][2] << " " <<
							all_norm[i][0] << " " << all_norm[i][1] << " " << all_norm[i][2] << "\n";
					}
				}

			}
			break;
		}
		case 1: {
			std::cout << "*** Convert external force to local coordinate ***\n";
			for (int i = frame0; i < frame1; i++) {

				int f = i * skipFactor;
				HermiteCurve curve;
				int seg_subdiv = 10;
				for (int y = 0; y < yarnNum; ++y) {

					//std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt"; //don't use upsampled centerline
					//const char* curvefile = tmp7.c_str();
					//std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";//don't use upsampled normals
					//const char* normfile = tmp8.c_str();
					std::string tmp9 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(y) + "_world.txt";
					const char* physical_world = tmp9.c_str();
					std::string tmp10 = "input/" + dataset + "/physical_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					const char* physical_local = tmp10.c_str();

					std::string tmp1 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* normfile_us = tmp1.c_str();
					std::string tmp2 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* curvefile_us = tmp2.c_str();

					std::string tmp5 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* twistfile = tmp5.c_str();

					std::ifstream fin2(curvefile_us);
					assert(fin2.is_open() && "curvefile_us file wasn't found!\n");
					//std::ifstream fin3(curvefile);
					//assert(fin3.is_open() && "curvefile file wasn't found!\n");
					std::ifstream fin4(physical_world);
					assert(fin4.is_open() && "physical_world file wasn't found!\n");

					curve.init_norm(curvefile_us, normfile_us, seg_subdiv); /******* normals are already twisted *******/
					std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;
					curve.assign(all_pts, all_tang, all_norm);
					//curve.assign_twist(twistfile, all_pts, all_tang, all_norm, 1);
					assert(all_pts.size() == yarn.getStepNum());			
					//curve.init_principleNormal(curvefile, normfile, seg_subdiv);


					//std::ifstream fin5(normfile);
					//assert(fin5.is_open() && "normfile file wasn't found!\n");  

					std::ifstream fin(physical_world);
					std::ofstream fout(physical_local);

					float S00, S01, S02, S10, S11, S12, S20, S21, S22;
					float A0, A1, A2;
					float B0, B1, B2;
					const int vrtx_num = yarn.getStepNum();
					
					for (int v = 0; v < vrtx_num; ++v) {
						fin >> S00 >> S01 >> S02 >> S10 >> S11 >> S12 >> S20 >> S21 >> S22
							//>> A0 >> A1 >> A2
							//>> B0 >> B1 >> B2
							;
						
						/*
						const double curveLength = curve.totalLength();
						float len = curveLength * (static_cast<double>(v) / static_cast<double>(vrtx_num - 1));
						const double t = curve.arcLengthInvApprox(len);
						Eigen::Vector3d ex, ey, ez;
						curve.getRotatedFrame(t, ex, ey, ez); */
						

						Eigen::Vector3d ez = all_tang[v];
						Eigen::Vector3d ey = all_norm[v];
						Eigen::Vector3d ex = ez.cross(ey);

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

						/// uncomment for when having the internal forces
						//Eigen::Vector3f localA, localB, worldA, worldB;
						//worldA << A0, A1, A2;
						//worldB << B0, B1, B2;
						//localA = M*worldA;
						//localB = M*worldB;
						//fout << localA(0) << " " << localA(1) << " " << localA(2) << " " <<
						//localB(0) << " " << localB(1) << " " << localB(2) ;

						fout << std::endl;
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
					//std::string tmp4 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
					//const char* curvefile = tmp4.c_str();
					//std::string tmp5 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
					//const char* normfile = tmp5.c_str();

					std::string tmp6 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* normfile_us = tmp6.c_str();
					std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* curvefile_us = tmp7.c_str();
					
					std::string tmp5 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* twistfile = tmp5.c_str();

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
					//const int upsample = 2;
					extractCompress_seg(configfile, yarnfile1, yarnfile2, "noNeed.txt", compress_S,
						curvefile_us, normfile_us, twistfile, 2, yarn.getPlyNum(), vrtx_num);
					/*************************************************/
					std::string tmp8 = "output/" + dataset + "/genYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					const char* outfile = tmp8.c_str();
					//// Procedural step
					yarn.simulate_ply();
					yarn.write_plys("test_ply.txt");
					const int K = yarn.getPlyNum();
					yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
					yarn.build("test_fly.txt", K);

					////pipeline 2:
					////yarn.compress_yarn3D(deformGrad, compress_S);

					yarn.compress_yarn_A(compress_S);
					yarn.curve_yarn(curvefile_us, normfile_us);
					yarn.write_yarn(outfile);
					
					///////*************************************************/
					//std::string tmp7 = "output/" + dataset + "/genYarn_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					//const char* outfile_wo = tmp7.c_str();
					//yarn.simulate_ply();
					//yarn.write_plys("test_ply.txt");
					//yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
					//yarn.build("test_fly.txt", K);
					//yarn.curve_yarn(curvefile, normfile);
					//yarn.write_yarn(outfile_wo);
				}
			}
			break;
		}
		case 2: {
			std::cout << "*** Training phase ***\n";

			for (int i = frame0; i < frame1; i++) {
				int f = i * skipFactor;
				for (int y = 0; y < yarnNum; ++y) {

					std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* curvefile_us = tmp1.c_str();
					std::string tmp2 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* normfile_us = tmp2.c_str();

					std::string tmp6 = "input/" + dataset + "/NN/testY_NN_full_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					const char* compress_S = tmp6.c_str();
					//std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
					//const char* curvefile = tmp7.c_str();
					//std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
					//const char* normfile = tmp8.c_str();

					//std::string tmp9 = "input/" + dataset + "/deformGrad_" + std::to_string(f) + "_trans.txt";
					//const char* deformGrad = tmp9.c_str();
					std::cout << compress_S << std::endl;
					std::ifstream fin2(compress_S);
					assert(fin2.is_open() && "compress_S_NN file wasn't found!\n");
					std::ifstream fin3(curvefile_us);
					assert(fin3.is_open() && "curvefile file wasn't found!\n");
					std::ifstream fin4(normfile_us);
					assert(fin4.is_open() && "normfile file wasn't found!\n");

					///*******  write the yarn ******/
					std::string tmp3 = "output/" + dataset + "/genYarn_NN_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					const char* outfile = tmp3.c_str();
					// Procedural step
					yarn.yarn_simulate();
					//pipeline 2:
					//yarn.compress_yarn3D(deformGrad, compress_S);
					yarn.compress_yarn_A(compress_S);
					yarn.curve_yarn(curvefile_us, normfile_us);
					yarn.write_yarn(outfile);
					//std::cout << outfile << std::endl;

					//*************** with fly-aways
					//std::string tmp4 = "output/" + dataset + "/genYarn_NN_flys_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					//const char* outfile_fly = tmp4.c_str();
					//// Procedural step
					//yarn.simulate_ply();
					//yarn.write_plys("test_ply.txt");
					//const int K = yarn.getPlyNum();
					//yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
					//yarn.build("test_fly.txt", K);
					//yarn.compress_yarn_A(compress_S);
					//yarn.curve_yarn(curvefile, normfile);
					//yarn.write_yarn(outfile_fly);
					//std::cout << outfile_fly << std::endl;

					/////*******  Validate NN by L2-norm ******/
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
			std::cout << "***Mapping yarn to a curve (generate fibers for frame0) *** \n";
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
				int seg_subdiv = 10;
				for (int y = 0; y < yarnNum; ++y) {

					//std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt"; //don't use upsampled centerline
					//const char* curvefile = tmp7.c_str();
					//std::string tmp8 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";//don't use upsampled normals
					//const char* normfile = tmp8.c_str();
					std::string tmp9 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(y) + "_world.txt";
					const char* physical_world = tmp9.c_str();
					std::string tmp10 = "input/" + dataset + "/physical_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					const char* physical_local = tmp10.c_str();

					std::string tmp1 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* normfile_us = tmp1.c_str();
					std::string tmp2 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* curvefile_us = tmp2.c_str();

					std::string tmp5 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* twistfile = tmp5.c_str();

					std::ifstream fin2(curvefile_us);
					assert(fin2.is_open() && "curvefile_us file wasn't found!\n");
					//std::ifstream fin3(curvefile);
					//assert(fin3.is_open() && "curvefile file wasn't found!\n");
					std::ifstream fin4(physical_world);
					assert(fin4.is_open() && "physical_world file wasn't found!\n");

					curve.init_norm(curvefile_us, normfile_us, seg_subdiv); /******* normals are already twisted *******/
					std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;
					curve.assign(all_pts, all_tang, all_norm);
					//curve.assign_twist(twistfile, all_pts, all_tang, all_norm, 1);
					assert(all_pts.size() == yarn.getStepNum());
					//curve.init_principleNormal(curvefile, normfile, seg_subdiv);


					//std::ifstream fin5(normfile);
					//assert(fin5.is_open() && "normfile file wasn't found!\n");  

					std::ifstream fin(physical_world);
					std::ofstream fout(physical_local);

					float S00, S01, S02, S10, S11, S12, S20, S21, S22;
					float A0, A1, A2;
					float B0, B1, B2;
					const int vrtx_num = yarn.getStepNum();

					for (int v = 0; v < vrtx_num; ++v) {
						fin >> S00 >> S01 >> S02 >> S10 >> S11 >> S12 >> S20 >> S21 >> S22
							//>> A0 >> A1 >> A2
							//>> B0 >> B1 >> B2
							;

						/*
						const double curveLength = curve.totalLength();
						float len = curveLength * (static_cast<double>(v) / static_cast<double>(vrtx_num - 1));
						const double t = curve.arcLengthInvApprox(len);
						Eigen::Vector3d ex, ey, ez;
						curve.getRotatedFrame(t, ex, ey, ez); */


						Eigen::Vector3d ez = all_tang[v];
						Eigen::Vector3d ey = all_norm[v];
						Eigen::Vector3d ex = ez.cross(ey);

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

						/// uncomment for when having the internal forces
						//Eigen::Vector3f localA, localB, worldA, worldB;
						//worldA << A0, A1, A2;
						//worldB << B0, B1, B2;
						//localA = M*worldA;
						//localB = M*worldB;
						//fout << localA(0) << " " << localA(1) << " " << localA(2) << " " <<
						//localB(0) << " " << localB(1) << " " << localB(2) ;

						fout << std::endl;
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
			skipFactor = 500;

			/**************** RUN ALL ****************/
			frame0 = 13000 / skipFactor ;		

			dataset = "spacing0.5x";
			frame1 = 14000 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.0x";
			frame1 = 14500 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.5x";
			frame1 = 15000 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			std::cout << " ********************** \n *********************** \n" << dataset << std::endl;

			dataset = "spacing0.5x_00011";
			frame1 = 16000 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing0.5x_10100";
			frame1 = 15000 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing0.5x_11110";
			frame1 = 15000 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1,yarnNum, dataset);
			std::cout << " ********************** \n *********************** \n" << dataset << std::endl;

			dataset = "spacing1.0x_00011";
			frame1 = 17000 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.0x_10100";
			frame1 = 15500 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.0x_11110";
			frame1 = 16000 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			std::cout << " ********************** \n *********************** \n" << dataset << std::endl;

			dataset = "spacing1.5x_00011";
			frame1 = 17500 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.5x_10100";
			frame1 = 16000 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.5x_11110";
			frame1 = 16500 / skipFactor + 1;
			phase1(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			std::cout << " ********************** \n *********************** \n" << dataset << std::endl;

			/********************************/
			break;
		}
		case 9: {
			/**************** RUN ALL ****************/
			frame0 = 8000 / skipFactor;

			dataset = "spacing0.5x";
			frame1 = 14000 / skipFactor + 1;
			//phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.0x";
			frame1 = 14500 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.5x";
			frame1 = 15000 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);

			dataset = "spacing0.5x_00011";
			frame1 = 16000 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing0.5x_10100";
			frame1 = 15000 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing0.5x_11110";
			frame1 = 15000 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.0x_00011";
			frame1 = 17000 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.0x_10100";
			frame1 = 15500 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.0x_11110";
			frame1 = 16000 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.5x_00011";
			frame1 = 17500 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.5x_10100";
			frame1 = 16000 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "spacing1.5x_11110";
			frame1 = 16500 / skipFactor + 1;
			phase2(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			/********************************/
			break;
		}
		case 10: {
			std::cout << "*** build traning data for one specific frame ***\n";

			std::string tmp0 = "input/" + dataset + "/centerYarn_" + std::to_string(15000) + "_" + std::to_string(0) + "_ds.txt";
			const char* curvefile_ds = tmp0.c_str();
			std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(15000) + "_" + std::to_string(0) + "_us.txt";
			const char* curvefile = tmp1.c_str();
			std::string tmp2 = "input/" + dataset + "/normYarn_" + std::to_string(15000) + "_" + std::to_string(0) + "_ds.txt";
			const char* normfile_ds = tmp2.c_str();
			std::string tmp3 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(15000) + "_" + std::to_string(0) + "_world.txt";
			const char* physical_world = tmp3.c_str();
			std::string tmp4 = "input/" + dataset + "/NN/trainX_" + std::to_string(15000) + "_" + std::to_string(0) + ".txt";
			const char* physical_local_window = tmp4.c_str();
			std::string tmp5 = "input/" + dataset + "/matrix_S_" + std::to_string(15000) + "_" + std::to_string(0) + ".txt";
			const char* compress_S = tmp5.c_str();
			std::string tmp6 = "input/" + dataset + "/NN/trainY_" + std::to_string(15000) + "_" + std::to_string(0) + ".txt";
			const char* compress_S_window = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/NN/angles_" + std::to_string(15000) + "_" + std::to_string(0) + ".txt";
			const char* angles = tmp7.c_str();

			

			std::ifstream fin1(curvefile_ds);
			assert(fin1.is_open() && "curvefile_ds file wasn't found!\n");
			std::ifstream fin3(physical_world);
			assert(fin3.is_open() && "physical_world file wasn't found!\n");
			std::ifstream fin4(compress_S);
			assert(fin4.is_open() && "compress_S file wasn't found!\n");

			std::ifstream fin_dg(physical_world);
			std::ifstream fin_S(compress_S);
			std::ofstream fout_dg(physical_local_window);
			std::ofstream fout_S(compress_S_window);
			std::ofstream fout_cntr(curvefile);
			std::ofstream fout_angle(angles);

			const int window_size = 9;
			int seg_subdiv = 10;
			HermiteCurve fullCurve;
			fullCurve.init(curvefile_ds, normfile_ds, seg_subdiv);
			const int vrtx_num = yarn.getStepNum();
			// write up-sampled centerline so later we can crop a segment out of it
			const double fullCurveLength = fullCurve.totalLength();
			std::vector<Eigen::Matrix3f> all_dg;
			std::vector<Eigen::Matrix2f> all_S;
			std::vector<Eigen::Vector3d> all_pnt;
			std::vector<Eigen::Vector3d> all_n;
			fout_cntr << vrtx_num << "\n";
			for (int v = 0; v < vrtx_num; ++v) {
				float fullLen = fullCurveLength * ( static_cast<double>(v) / static_cast<double>(vrtx_num - 1) );
				const double t_fll = fullCurve.arcLengthInvApprox(fullLen);
				Eigen::Vector3d pnt = fullCurve.eval(t_fll);
				fout_cntr << pnt[0] << " " << pnt[1] << " " << pnt[2] << "\n";

				//store normals for all points
				Eigen::Vector3d n = fullCurve.evalNormal(t_fll);
				all_n.push_back(n);
				/* Note that normals don't exactly match with up-sampled curve because adding new vertices to curve changes its curvature a bit */

				//store dg for all points
				float dg00, dg01, dg02, dg10, dg11, dg12, dg20, dg21, dg22;
				fin_dg >> dg00 >> dg01 >> dg02 >> dg10 >> dg11 >> dg12 >> dg20 >> dg21 >> dg22;
				Eigen::Matrix3f world_dg;
				world_dg << dg00, dg01, dg02,
					dg10, dg11, dg12,
					dg20, dg21, dg22;
				all_dg.push_back(world_dg);

				//store S-matrix for all points
				float S00, S01, S10, S11;
				fin_S >> S00 >> S01 >> S10 >> S11;
				Eigen::Matrix2f S;
				S << S00, S01, S10, S11;
				all_S.push_back(S);

			}
			fout_cntr.close();

			const float trimPercent = 0.10;
			const int ignorPlanes = trimPercent * vrtx_num; // crop the first and last 10% of the yarn
			for (int w = ignorPlanes; w < (vrtx_num - window_size + 1) - ignorPlanes; w++) {
				//define a curve segment 
				const int start = w;
				const int end = w + (window_size - 1);
				HermiteCurve curve;
				curve.init_seg(curvefile, start, end, seg_subdiv);				
				const double curveLength = curve.totalLength();

				for (int v = 0; v < window_size; ++v) {

					const double curveLength = curve.totalLength();
					float len = curveLength * (static_cast<double>(v) / static_cast<double>(window_size - 1));
					const double t = curve.arcLengthInvApprox(len);

					Eigen::Vector3d ex, ey, ez;
					curve.getRotatedFrame(t, ex, ey, ez);

					/** local to world **/
					Eigen::Matrix3f local_dg, M;
					M << ex[0], ex[1], ex[2],
						ey[0], ey[1], ey[2],
						ez[0], ez[1], ez[2];
					const int indx = w + v;
					local_dg = M*all_dg[indx]*M.transpose();

					//write converted parameters
					fout_dg << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
						local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
						local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
				}
				fout_dg << std::endl;

				const int v_full = ceil((start + end) / 2.0); //index for the full curve			
				Eigen::Vector3d n_full = all_n[v_full];

				const int v = ceil((end - start) / 2.0);
				float len = curveLength * (static_cast<double>(v) / static_cast<double>(window_size - 1));
				const double t = curve.arcLengthInvApprox(len);
				Eigen::Vector3d n = curve.evalNormal(t);

				//std::cout << " v: " << v_full << " " << v + start << std::endl;
				assert(v_full == v + start && "index for full yarn must be equal to index segment added with starting vertex");

				// rotate the shape-matching matrix to align the new normal
				const float angle = acos(n_full.dot(n));
				//std::cout << " dot product " << n_full.dot(n) << " angle " << angle << std::endl;
				Eigen::Matrix2f R, S, S_rot;
				R << cos(angle), -sin(angle),
					sin(angle), cos(angle);

				S_rot = R*all_S[v_full]*R.transpose();
				fout_S << S_rot(0, 0) << " " << S_rot(0, 1) << " " << S_rot(1, 0) << " " << S_rot(1, 1) << "\n";
				//std::cout << S << std::endl << S_rot << std::endl << std::endl;
				fout_angle << angle << std::endl;
			}

			fout_dg.close();
			fout_S.close();
			fout_angle.close();

			break;
		}
		case 11: {
			std::cout << "*** build training data for all frames ***\n";
			const float trimPercent = 0.15;
			const int window_size = 40;

			std::string tmp7 = "input/" + dataset + "/NN/trainX_all.txt";
			const char* all_trainX = tmp7.c_str();
			std::string tmp8 = "input/" + dataset + "/NN/trainY_all.txt";
			const char* all_trainY = tmp8.c_str();
			std::ofstream fout_trainX_all(all_trainX);
			std::ofstream fout_trainY_all(all_trainY);

			for (int i = frame0; i < frame1; i++) {

				int f = i * skipFactor;
				for (int y = 0; y < yarnNum; ++y) {

					std::string tmp0 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
					const char* curvefile_ds = tmp0.c_str();
					std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
					const char* curvefile = tmp1.c_str();
					std::string tmp2 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_ds.txt";
					const char* normfile_ds = tmp2.c_str();
					std::string tmp3 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(y) + "_world.txt";
					const char* physical_world = tmp3.c_str();
					std::string tmp4 = "input/" + dataset + "/NN/trainX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					const char* physical_local_window = tmp4.c_str();
					std::string tmp5 = "input/" + dataset + "/matrix_S_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					const char* compress_S = tmp5.c_str();
					std::string tmp6 = "input/" + dataset + "/NN/trainY_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					const char* compress_S_window = tmp6.c_str();
					std::string tmp7 = "input/" + dataset + "/NN/angles_" + std::to_string(f) + "_" + std::to_string(0) + ".txt";
					const char* angles = tmp7.c_str();

					std::string tmp8 = "input/" + dataset + "/NN/testX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
					const char* physical_local_window_test = tmp8.c_str();

					buildTraining(curvefile_ds, normfile_ds, physical_world, compress_S, yarn, trimPercent, window_size, curvefile, angles, 
						physical_local_window, compress_S_window, physical_local_window_test, fout_trainX_all, fout_trainY_all, 1);
				}
			}

			fout_trainX_all.close();
			fout_trainY_all.close();

			break;
		}
		case 12: {

			int yarnNum = 1;
			int skipFactor = 100;
			int frame0 = 200 / skipFactor;
			float trimPercent = 0;
			int window_size = 50;
			int isTrain = 0;

			int frame1 = 200 / skipFactor + 1;
			std::string dataset = "spacing1.0x_00011_woven";
			//buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

			///////////////////////
			trimPercent = 0;
			yarnNum = 1;
			skipFactor = 500;
			frame0 = 17000 / skipFactor;
			window_size = 50;
			isTrain = 1;

			//int frame1 = 16000 / skipFactor + 1;
			//std::string dataset = "spacing0.5x_00011";
			//buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

			//frame1 = 15000 / skipFactor + 1;
			//dataset = "spacing0.5x_10100";
			//buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

			//frame1 = 15000 / skipFactor + 1;
			//dataset = "spacing0.5x_11110";
			//buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

			frame1 = 17000 / skipFactor + 1;
			dataset = "spacing1.0x_00011";
			buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

			//frame1 = 15500 / skipFactor + 1;
			//dataset = "spacing1.0x_10100";
			//buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

			//frame1 = 16000 / skipFactor + 1;
			//dataset = "spacing1.0x_11110";
			//buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

		    //frame1 = 17500 / skipFactor + 1;
			//dataset = "spacing1.5x_00011";
			//buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

			//frame1 = 16000 / skipFactor + 1;
			//dataset = "spacing1.5x_10100";
			//buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

			//frame1 = 16500 / skipFactor + 1;
			//dataset = "spacing1.5x_11110";
			//buildTraning_all(yarn, skipFactor, frame0, frame1, yarnNum, dataset, window_size, trimPercent, isTrain);

			break;
		}
		case 13: {
			std::cout << "*** Build training data - new *** \n";

			std::string tmp0 = "input/" + dataset + "/NN/trainX_all.txt";
			const char* all_trainX = tmp0.c_str();
			std::string tmp1 = "input/" + dataset + "/NN/trainY_all.txt";
			const char* all_trainY = tmp1.c_str();
			std::ofstream fout_trainX_all(all_trainX);
			std::ofstream fout_trainY_all(all_trainY);

			int f = 200;

			//std::string tmp0 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(0) + "_ds.txt";
			//const char* curvefile_ds = tmp0.c_str();
			std::string tmp2 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(0) + "_ds.txt";
			const char* normfile_ds = tmp2.c_str();
			std::string tmp3 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(0) + "_world.txt";
			const char* physical_world = tmp3.c_str();
			std::string tmp4 = "input/" + dataset + "/NN/trainX_" + std::to_string(f) + "_" + std::to_string(0) + ".txt";
			const char* physical_local_seg = tmp4.c_str();
			std::string tmp5 = "input/" + dataset + "/matrix_S_" + std::to_string(f) + "_" + std::to_string(0) + ".txt";
			const char* compress_S = tmp5.c_str();
			std::string tmp6 = "input/" + dataset + "/NN/trainY_" + std::to_string(f) + "_" + std::to_string(0) + ".txt";
			const char* compress_S_seg = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/NN/angles_" + std::to_string(f) + "_" + std::to_string(0) + ".txt";
			const char* angles = tmp7.c_str();

			std::string tmp8 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(0) + "_us.txt";
			const char* twistfile = tmp8.c_str();
			std::string tmp9 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(0) + "_us.txt";
			const char* curvefile_us = tmp9.c_str();
			std::string tmp10 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(0) + "_us.txt";
			const char* normfile_us = tmp10.c_str();


			//std::ifstream fin1(curvefile_ds);
			//assert(fin1.is_open() && "curvefile_ds file wasn't found!\n");
			std::ifstream fin3(physical_world);
			assert(fin3.is_open() && "physical_world file wasn't found!\n");

			const int seg_subdiv = 10;
			const int window_size = 50;
			const int isTrain = 0;
			const float trimPercent = 0.0;
			const int vrtx_num = yarn.getStepNum();

			std::ofstream fout_trainX(physical_local_seg);
			std::ofstream fout_trainY(compress_S_seg);
			std::ofstream fout_angle(angles);

			/* yarn-level */
			std::vector<Eigen::Matrix3f> all_dg;
			assign_dg(physical_world, all_dg);
			std::vector<Eigen::Matrix2f> all_S;
			if (isTrain) {
				assign_S(compress_S, all_S);
			}
			std::vector<float> twists;
			assign_twist(twistfile, twists);

			HermiteCurve curve;
			curve.init_norm(curvefile_us, normfile_us, seg_subdiv);
			std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;
			curve.assign(all_pts, all_tang, all_norm);
			//curve.assign_twist(twists, all_pts, all_tang, all_norm, 2);

			/* window-level */
			const int ignorPlanes = trimPercent * vrtx_num; // crop the first and last #% of the yarn
			for (int w = ignorPlanes; w < (vrtx_num - window_size + 1) - ignorPlanes; w++) {
				//std::cout << w << std::endl;

				//define a curve segment 
				const int start = w;
				const int end = w + (window_size - 1);

				HermiteCurve segment;
				segment.init_seg(curvefile_us, start, end, seg_subdiv);
				std::vector<Eigen::Vector3d> all_pts_seg, all_tang_seg, all_norm_seg;
				//segment.assign(all_pts_seg, all_tang_seg, all_norm_seg);
				std::vector<float>::const_iterator first = twists.begin() + start;
				std::vector<float>::const_iterator last = twists.begin() + end + 1;
				std::vector<float> twists_seg(first, last);
				segment.assign_twist(twists_seg, all_pts_seg, all_tang_seg, all_norm_seg, 1);

				std::vector<Eigen::Matrix3f> all_dg_seg, all_local_dg_seg;
				for (int d = start; d <= end; d++) all_dg_seg.push_back(all_dg[d]);
				transfer_dg_2local(all_tang_seg, all_norm_seg, all_dg_seg, all_local_dg_seg); //later add augmentation



				for (int d = 0; d < window_size; d++) {
					Eigen::Matrix3f local_dg = all_local_dg_seg[d];
					//std::cout << all_local_dg[d] << std::endl << std::endl;
					fout_trainX << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
						local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
						local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
					fout_trainX_all << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
						local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
						local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
				}
				fout_trainX << "\n";
				fout_trainX_all << "\n";

				const int v_yarn = ceil((start + end) / 2.0);
				const int v_seg = ceil((end - start) / 2.0);
				assert(v_yarn == v_seg + start && "index for full yarn must be equal to index segment added with starting vertex");

				Eigen::Vector3d norm1 = all_norm[v_yarn];
				Eigen::Vector3d norm2 = all_norm_seg[v_seg];
				Eigen::Vector3d tang = all_tang[v_yarn];
				assert((all_tang[v_yarn] - all_tang_seg[v_seg]).norm() < eps && "tangents for both frames must be similar!\n");
				float angle = get_angle(norm1, norm2, tang);
				fout_angle << angle << "\n";

				if (isTrain) {
					Eigen::Matrix2f S_local;
					rotate_S_2local(all_S[v_yarn], S_local, angle);
					fout_trainY << S_local(0, 0) << " " << S_local(0, 1) << " " << S_local(1, 0) << " " << S_local(1, 1) << "\n";
					fout_trainY_all << S_local(0, 0) << " " << S_local(0, 1) << " " << S_local(1, 0) << " " << S_local(1, 1) << "\n";
				}
			}
			fout_trainX_all.close();
			fout_trainY_all.close();
			fout_trainX.close();
			fout_trainY.close();
			fout_angle.close();


			//Eigen::Vector3d cross = all_norm[61].cross(all_norm_seg[1]);
			//std::cout << "seg --- \n " << all_norm_seg[2].dot(all_tg_seg[2]) << std::endl << all_norm[53].dot(all_tg[53]) << std::endl <<  std::endl;
			//std::cout << "seg --- \n " << cross  << std::endl << all_tg_seg[1] << std::endl << all_tg[61] << std::endl;
			//std::cout << "full ********** \n " << all_tang_seg[v_seg] << std::endl << all_norm[2] << std::endl << all_tang[v_yarn] << std::endl;


			//std::ofstream fout_TNB("../data/TNB.txt");
			//for (int i = 50; i < 100; i++) {
			//	fout_TNB << all_pts[i][0] << " " << all_pts[i][1] << " " << all_pts[i][2] << " " <<
			//		all_tg[i][0] << " " << all_tg[i][1] << " " << all_tg[i][2] << " " <<
			//		all_norm[i][0] << " " << all_norm[i][1] << " " << all_norm[i][2] << "\n";
			//}
			break;
		}
		case 100: {
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