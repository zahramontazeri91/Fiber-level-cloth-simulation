#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"
#include "curveFitting.h"

#include <stdio.h>
#include <cmath>


void writeParameters(std::vector<Matrix_S> &all_mat_S, std::vector<float> &all_theta_R, const char* compress_R, const char* compress_S)
{

	FILE *foutR;
	if (fopen_s(&foutR, compress_R, "wt") == 0) {
		fprintf_s(foutR, "%d \n", all_theta_R.size());
		for (int i = 0; i < all_theta_R.size(); ++i) {
			fprintf_s(foutR, "%.6f \n", all_theta_R[i]);
		}
		fclose(foutR);
	}

	FILE *foutS;
	// write S-matrix for each segment not vertex 
	if (fopen_s(&foutS, compress_S, "wt") == 0) {
		//fprintf_s(foutS, "%d \n", all_mat_S.size());
		for (int i = 0; i < all_mat_S.size() - 1; ++i) {
			fprintf_s(foutS, "%.6f %.6f %.6f \n", (all_mat_S[i].S00 + all_mat_S[i + 1].S00) / 2.f,
				(all_mat_S[i].S11 + all_mat_S[i + 1].S11) / 2.f,
				(all_mat_S[i].S01 + all_mat_S[i + 1].S01) / 2.f);
		}
		int i = all_mat_S.size() - 1;
		fprintf_s(foutS, "%.6f %.6f %.6f \n", all_mat_S[i].S00 / 2.f,
			all_mat_S[i].S11 / 2.f,
			all_mat_S[i].S01 / 2.f);
		fclose(foutS);
	}
	std::cout << "Compression parameters are successfully written to the file!\n";
}
void decomposeS(const Matrix_S &mat_S, Ellipse &ellipse) {
	Eigen::Matrix2f S;
	S << mat_S.S00, mat_S.S01, mat_S.S01, mat_S.S11;
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> eigensolver(S);
	Eigen::Matrix2f sigma = eigensolver.eigenvalues().asDiagonal();
	Eigen::Matrix2f V = eigensolver.eigenvectors();

	// S = V*sigma*V_T
	ellipse.longR = std::max(sigma(0, 0), sigma(1, 1));
	ellipse.shortR = std::min(sigma(0, 0), sigma(1, 1));
	//Make V reflection free (if det(V) is negative)
	Eigen::Matrix2f m;
	m << 1, 0, 0, -1;
	if (V.determinant() < 0)
		V = V*m;
	ellipse.angle = atan2(V(1, 0), V(0, 0));
	ellipse.angle = ellipse.angle < 0 ? ellipse.angle + 2.f*pi : ellipse.angle;
}


void extractCompress_seg(const char* configfile, const char* yarnfile1, const char* yarnfile2, const char* deformGrad, const char* compress_S,
	const char* curveFile, const char* normFile, const char* global_rot, const int ply_num, const int vrtx_num)
{
	const int n = vrtx_num;

	//pipeline 1:
	std::vector<yarnIntersect2D> pnts_ref;
	CrossSection cs1(yarnfile1, configfile, pnts_ref);

	//Fiber::Yarn yarn_tmp;
	//yarn_tmp.yarnCenter(yarnfile2, curveFile);
	std::vector<yarnIntersect2D> pnts_trans;
	CrossSection cs2(yarnfile2, curveFile, normFile, ply_num, n, 100, pnts_trans, true);
	//CrossSection cs2(upsample, yarnfile2, curveFile, normFile, twistFile, ply_num, n, 100, pnts_trans, true);  //changed this to true 

	//write global rotations for phase-matching purpose
	std::ofstream phase_fout(global_rot);
	std::vector<Eigen::Matrix2f> all_A;
	cs2.yarnShapeMatches_A(pnts_trans, pnts_ref, all_A, phase_fout);
	phase_fout.close();

	FILE *foutS;
	// write S-matrix for each segment not vertex 
	if (fopen_s(&foutS, compress_S, "wt") == 0) {
		//fprintf_s(foutS, "%d \n", all_mat_S.size());
		for (int i = 0; i < all_A.size(); ++i) {
			fprintf_s(foutS, "%.6f %.6f %.6f %.6f \n", all_A[i](0, 0), all_A[i](0, 1), all_A[i](1, 0), all_A[i](1, 1));
		}
		fclose(foutS);
	}

	//constant fitting
	//std::cout << "Compression parameters are successfully written to the file!\n";

	//for debug: visualization
	//const char* L2File = "../data/L2.txt";
	//const char* refFile = "../data/allCrossSection2D_ref.txt";
	//const char* deformedRefFile = "../data/allCrossSection2D_deformedRef.txt";
	//const char* deformedFile = "../data/allCrossSection2D_deformed.txt";
	//const float trimPercent = 0.2;
	//plotIntersections(pnts_ref, refFile, trimPercent);
	//std::vector<yarnIntersect2D> ref_deformed;
	//(pnts_ref, ref_deformed, all_A);
	//plotIntersections(ref_deformed, deformedRefFile, trimPercent);
	//plotIntersections(pnts_trans, deformedFile, trimPercent);

	//std::vector<float> L2;
	//L2norm(ref_deformed, pnts_trans, L2, L2File, trimPercent); //note that these have same size
}

float nextTheta(float theta0, float theta1) {
	const int k = floor((theta0 - theta1)*0.5 / pi);
	float ans = theta1 + (k - 1)*2.0*pi;
	float best = abs(theta0 - ans);
	for (int i = 0; i < 2; ++i) {
		float cur = theta1 + (k + i)*2.0*pi;
		float val = abs(theta0 - cur);
		if (best > val) {
			best = val;
			ans = cur;
		}
	}
	return ans;
}

void nonPeriodicTheta(const std::vector<float> &theta, std::vector<float> &theta_new) {
	const int m = theta.size();
	theta_new.resize(m);
	theta_new = theta;
	for (int i = 0; i < m; ++i) {
		theta_new[i] = nextTheta(theta_new[i - 1], theta_new[i]);
	}
	float k = floor(0.25*(theta_new[0] + theta_new[m - 1]) / pi + 0.5)*2.0*pi;
	for (int i = 0; i < m; ++i) {
		theta_new[i] -= k;
	}
}

void L2norm(const std::vector<yarnIntersect2D> &its_deform, const std::vector<yarnIntersect2D> &its_trans, std::vector<float> &L2, const char* filename, const float trimPercent) {

	if (its_deform.size() != its_trans.size())
		std::cout << its_deform.size() << " " << its_trans.size() << std::endl;
	assert(its_deform.size() == its_trans.size());
	FILE *fout;

	const int ignorPlanes = trimPercent * its_deform.size(); // crop the first and last 10% of the yarn
	if (fopen_s(&fout, filename, "wt") == 0) {
		fprintf_s(fout, "%d \n", its_deform.size());
		for (int i = ignorPlanes; i < its_deform.size() - ignorPlanes; ++i) { //number of planes
			float e = 0.f;
			for (int p = 0; p < its_deform[i].size(); ++p) { //number of plys
				for (int j = 0; j < its_deform[i][p].size(); ++j) { //number of intersections
					e += square_norm(its_deform[i][p][j] - its_trans[i][p][j]);
				}
			}
			L2.push_back(e);
			fprintf_s(fout, "%.6f \n", e);
		}
		fclose(fout);
	}
}

void plotIntersections(const std::vector<yarnIntersect2D> &its, const char* filename, const float trimPercent) {
	FILE *fout;
	// write the plycenters
	if (fopen_s(&fout, filename, "wt") == 0) {
		const int ignorPlanes = trimPercent * its.size(); // crop the first and last 10% of the yarn

		fprintf_s(fout, "plane_num: %d \n", its.size() - 2 * ignorPlanes);
		fprintf_s(fout, "ply_num: %d \n", its[0].size());
		fprintf_s(fout, "\n");

		for (int i = ignorPlanes; i < its.size() - ignorPlanes; ++i) { //number of planes
			fprintf_s(fout, "index_plane : %d \n", i - ignorPlanes);
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

void deformRef(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Eigen::Matrix2f> &all_A) {

	its_deformed.resize(its.size());
	Eigen::Matrix2f transf;
	for (int i = 0; i < its_deformed.size(); ++i) { //number of planes
		transf << all_A[i](0, 0), all_A[i](0, 1), all_A[i](1, 0), all_A[i](1, 1);

		its_deformed[i].resize(its[i].size());
		for (int p = 0; p < its[i].size(); ++p) { //number of plys
			its_deformed[i][p].resize(its[i][p].size());
			for (int j = 0; j < its[i][p].size(); ++j) { //number of intersections
				Eigen::VectorXf ref(2, 1);
				ref << its[i][p][j].x, its[i][p][j].y;
				Eigen::VectorXf deformedRef = transf*ref;
				its_deformed[i][p][j] = vec2f(deformedRef(0, 0), deformedRef(1, 0));
			}
		}
	}
}

void deformRef(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Matrix_S> &all_mat_S, const std::vector<float> &all_theta_R) {

	its_deformed.resize(its.size());
	Eigen::Matrix2f R, S, V, sigma, transf;
	for (int i = 0; i < its_deformed.size(); ++i) { //number of planes
		float S00 = all_mat_S[i].S00;
		float S01 = all_mat_S[i].S01;
		float S11 = all_mat_S[i].S11;
		Eigen::Matrix2f S;
		S << S00, S01, S01, S11;
		R << cos(all_theta_R[i]), -sin(all_theta_R[i]), sin(all_theta_R[i]), cos(all_theta_R[i]);
		transf = R *S;

		its_deformed[i].resize(its[i].size());
		for (int p = 0; p < its[i].size(); ++p) { //number of plys
			its_deformed[i][p].resize(its[i][p].size());
			for (int j = 0; j < its[i][p].size(); ++j) { //number of intersections
				Eigen::VectorXf ref(2, 1);
				ref << its[i][p][j].x, its[i][p][j].y;
				Eigen::VectorXf deformedRef = transf*ref;
				its_deformed[i][p][j] = vec2f(deformedRef(0, 0), deformedRef(1, 0));
			}
		}
	}
}

void deformRef(const std::vector<yarnIntersect2D> &its, std::vector<yarnIntersect2D> &its_deformed,
	const std::vector<Ellipse> &ellipses, const std::vector<float> &all_theta_R) {

	its_deformed.resize(its.size());
	Eigen::Matrix2f R, S, V, sigma, transf;
	for (int i = 0; i < its_deformed.size(); ++i) { //number of planes
		sigma << ellipses[i].longR, 0, 0, ellipses[i].shortR;
		V << cos(ellipses[i].angle), -sin(ellipses[i].angle), sin(ellipses[i].angle), cos(ellipses[i].angle);
		S = V*sigma*V.transpose();
		R << cos(all_theta_R[i]), -sin(all_theta_R[i]), sin(all_theta_R[i]), cos(all_theta_R[i]);
		transf = R *S;

		its_deformed[i].resize(its[i].size());
		for (int p = 0; p < its[i].size(); ++p) { //number of plys
			its_deformed[i][p].resize(its[i][p].size());
			for (int j = 0; j < its[i][p].size(); ++j) { //number of intersections
				Eigen::VectorXf ref(2, 1);
				ref << its[i][p][j].x, its[i][p][j].y;
				Eigen::VectorXf deformedRef = transf*ref;
				its_deformed[i][p][j] = vec2f(deformedRef(0, 0), deformedRef(1, 0));
			}
		}
	}
}

void assign_dg(const char* physical_world, std::vector<Eigen::Matrix3f> &all_world_dg) {
	std::ifstream fin(physical_world);
	assert(fin.is_open() && "physical_world file wasn't found!\n");

	//store dg for all points
	while (1) {
		float dg00, dg01, dg02, dg10, dg11, dg12, dg20, dg21, dg22;
		fin >> dg00 >> dg01 >> dg02 >> dg10 >> dg11 >> dg12 >> dg20 >> dg21 >> dg22;
		if (fin.eof()) break;
		Eigen::Matrix3f world_dg;
		world_dg << dg00, dg01, dg02,
			dg10, dg11, dg12,
			dg20, dg21, dg22;
		all_world_dg.push_back(world_dg);
	}

	fin.close();
}

void assign_S(const char* compress_S, std::vector<Eigen::Matrix2f> &all_S) {
	std::ifstream fin(compress_S);
	assert(fin.is_open() && "compress_S file wasn't found!\n");

	//store dg for all points
	while (1) {
		//store S-matrix for all points
		float S00, S01, S10, S11;
		fin >> S00 >> S01 >> S10 >> S11;
		if (fin.eof()) break;
		Eigen::Matrix2f S;
		S << S00, S01, S10, S11;
		all_S.push_back(S);
	}

	fin.close();
}

void assign_twist(const char* twistFile, std::vector<float> &twists) {
	std::ifstream fin(twistFile);
	assert(fin.is_open() && "twist-file is not found!\n");

	int n = 0;
	fin >> n;
	for (int i = 0; i < n; i++) {
		float twist;
		fin >> twist;
		twists.push_back(twist);
	}
}

void transfer_dg_2local(const std::vector<Eigen::Vector3d> &all_tang, const std::vector<Eigen::Vector3d> &all_norm,
	const std::vector<Eigen::Matrix3f> &all_world_dg, std::vector<Eigen::Matrix3f> &all_local_dg, const int flip) {

	assert(all_tang.size() == all_world_dg.size() && all_norm.size() == all_world_dg.size());
	const int n = all_world_dg.size();
	Eigen::Vector3d ex, ey, ez;
	for (int i = 0; i < n; i++) {


		if (flip == 1) { //flip normal
			ez = all_tang[i];
			ey = -1.0*all_norm[i];
		}
		else if (flip == 2) { //flip tangent
			ez = -1.0*all_tang[i];
			ey = all_norm[i];
		}
		else if (flip == 3) { //flip tangent
			ez = -1.0*all_tang[i];
			ey = -1.0*all_norm[i];
		}
		else {
			ez = all_tang[i];
			ey = all_norm[i];
		}
		ex = ez.cross(ey);

		/** world to local **/
		Eigen::Matrix3f local_dg, M;
		M << ez[0], ez[1], ez[2],
			ey[0], ey[1], ey[2],
			ex[0], ex[1], ex[2];

		////local_dg = M*all_world_dg[i] * M.transpose();
		local_dg = M*all_world_dg[i];

		all_local_dg.push_back(local_dg);
	}
}

void rotate_S_2local(const Eigen::Matrix2f &S, Eigen::Matrix2f &S_local, const float &angle, const int flip) {

	Eigen::Matrix2f R, A;

	R << cos(angle), -sin(angle),
		sin(angle), cos(angle);
	S_local = R*S* R.transpose();

	if (flip == 1) //flip normal
		A << -1, 0, 0, -1;
	else if (flip == 2)  //flip tangent
		A << 1, 0, 0, -1;
	else if (flip == 3)  //flip tangent and normal
		A << -1, 0, 0, -1;
	else
		A << 1, 0, 0, 1;

	S_local = A*S_local*A.transpose();

}

float get_angle(Eigen::Vector3d &norm1, Eigen::Vector3d &norm2, Eigen::Vector3d &tang) {

	// rotate the shape-matching matrix to align the new normal
	float angle = acos(norm1.dot(norm2));
	if (angle != angle) //dot product (must be 1) but might be larger than 1 and so acos return nan 
		angle = norm1.dot(norm2) > 0 ? 0.0 : pi;

	Eigen::Vector3d cross = (norm1.cross(norm2)).normalized();
	if (cross.norm() < eps) return 0.0;

	float dist = (cross - tang).norm();
	float dist_ = (cross - (-1.0*tang)).norm();
	angle = dist < dist_ ? angle : -1.0*angle;

	return angle;
}


void loadSamples(const char* curvefile, std::vector<Eigen::Vector3f> &pnts) {
	std::ifstream fin(curvefile);
	if (!fin.is_open())
		std::cout << curvefile << std::endl;
	assert(fin.is_open() && "curvefile file wasn't found!\n");
	int pnts_num = 0;
	fin >> pnts_num;
	Eigen::Vector3f pnt;
	for (int i = 0; i < pnts_num; i++) {
		fin >> pnt[0] >> pnt[1] >> pnt[2];
		pnts.push_back(pnt);
	}
}

void findAABB(std::vector<Eigen::Vector3f> &pnts, float &min_x, float &min_y, float &min_z, float &max_x, float &max_y, float &max_z) {

	const int sz = pnts.size();
	max_x = std::numeric_limits<float>::min();
	max_y = max_x;
	max_z = max_x;
	min_x = std::numeric_limits<float>::max();
	min_y = min_x;
	min_z = min_x;
	for (int i = 0; i < sz; i++) {
		//minimum
		if (pnts[i][0] < min_x)
			min_x = pnts[i][0];
		if (pnts[i][1] < min_y)
			min_y = pnts[i][1];
		if (pnts[i][2] < min_z)
			min_z = pnts[i][2];
		// maximum
		if (pnts[i][0] > max_x)
			max_x = pnts[i][0];
		if (pnts[i][1] > max_y)
			max_y = pnts[i][1];
		if (pnts[i][2] > max_z)
			max_z = pnts[i][2];
	}

	//std::cout << min_x << " " << min_y << " " << min_z << " " << max_x << " " << max_y << " " << max_z;
}

void fillVolume(const std::vector<Eigen::Vector3f> &pnts, const float radius, const float minAABB[3], const float maxAABB[3], const int resol[3], std::vector<std::vector<std::vector<float>>> &vol) {

	//initialize vol
	vol.resize(resol[0]);
	for (int x = 0; x < resol[0]; x++) {
		vol[x].resize(resol[1]);
		for (int y = 0; y < resol[1]; y++) {
			vol[x][y].resize(resol[2]);
			for (int z = 0; z < resol[2]; z++) {
				vol[x][y].push_back(0.f);
			}
		}
	}

	const int sz = pnts.size();
	for (int i = 0; i < sz; i++) {
		const float len_x = maxAABB[0] - minAABB[0];
		const float len_y = maxAABB[1] - minAABB[1];
		const float len_z = maxAABB[2] - minAABB[2];

		int idx_x = ((pnts[i][0] - minAABB[0]) / len_x) * resol[0];
		int idx_y = ((pnts[i][1] - minAABB[1]) / len_y) * resol[1];
		int idx_z = ((pnts[i][2] - minAABB[2]) / len_z) * resol[2];

		if (idx_x == resol[0]) idx_x = idx_x - 1;
		if (idx_y == resol[1]) idx_y = idx_y - 1;
		if (idx_z == resol[2]) idx_z = idx_z - 1;

		vol[idx_x][idx_y][idx_z] = 1.f;

		// go d distance in all 6 directions and add neighbor voxels if needed
		float bottom_x = minAABB[0] + idx_x * (len_x / float(resol[0]));
		float top_x = bottom_x + (len_x / float(resol[0]));

		float bottom_y = minAABB[1] + idx_y * (len_y / float(resol[1]));
		float top_y = bottom_y + (len_y / float(resol[1]));

		float bottom_z = minAABB[2] + idx_z * (len_z / float(resol[2]));
		float top_z = bottom_z + (len_z / float(resol[2]));


		if ((top_x - pnts[i][0]) < radius)
			if (idx_x + 1 != resol[0])
				vol[idx_x + 1][idx_y][idx_z] = 1.f;
		if ((pnts[i][0] - bottom_x) < radius)
			if (idx_x - 1 >= 0)
				vol[idx_x - 1][idx_y][idx_z] = 1.f;

		if ((top_y - pnts[i][1]) < radius)
			if (idx_y + 1 != resol[1])
				vol[idx_x][idx_y + 1][idx_z] = 1.f;
		if ((pnts[i][1] - bottom_y) < radius)
			if (idx_y - 1 >= 0)
				vol[idx_x][idx_y - 1][idx_z] = 1.f;

		if ((top_z - pnts[i][2]) < radius)
			if (idx_z + 1 != resol[2])
				vol[idx_x][idx_y][idx_z + 1] = 1.f;
		if ((pnts[i][2] - bottom_z) < radius)
			if (idx_z - 1 >= 0)
				vol[idx_x][idx_y][idx_z - 1] = 1.f;

	}
}

void writeVol(std::string &dataset, const int frame, const int yarn0, const int yarn1, const int resol_x, const int resol_y, const int resol_z, const float radius, const char* volumeFile) {

	std::vector<Eigen::Vector3f> pnts;
	//for loop over yarns y
	for (int y = yarn0; y < yarn1; y++) {
		std::string tmp = "input/" + dataset + "/centerYarn_" + std::to_string(frame) + "_" + std::to_string(y) + "_us.txt";
		const char* curvefile_us = tmp.c_str();
		loadSamples(curvefile_us, pnts);
	}

	float min_x, min_y, min_z, max_x, max_y, max_z;
	findAABB(pnts, min_x, min_y, min_z, max_x, max_y, max_z);


	float minAABB[3], maxAABB[3];
	float *data;
	float *data_vol;
	//int resol[3];
	int N;


	minAABB[0] = min_x - radius, minAABB[1] = min_y - radius, minAABB[2] = min_z - radius;
	maxAABB[0] = max_x + radius, maxAABB[1] = max_y + radius, maxAABB[2] = max_z + radius;


	// Modify tile scale here
	float scale = 1.0;

	int resol[3];
	resol[0] = resol_x;
	resol[1] = resol_y;
	resol[2] = resol_z;

	N = resol[0] * resol[1] * resol[2];
	data = new float[N];
	data_vol = new float[N];


	std::vector<std::vector<std::vector<float> > > volume;
	fillVolume(pnts, radius, minAABB, maxAABB, resol, volume);
	// flatten the volume
	int i = 0;
	for (int z = 0; z < resol[2]; z++) {
		for (int y = 0; y < resol[1]; y++) {
			for (int x = 0; x < resol[0]; x++) {
				data_vol[i] = volume[x][y][z];
				i++;
			}
		}
	}

	//Modidy data here
	for (int i = 0; i < N; i++) {
		data[i] = data_vol[i];
		//data[i] = 0.1;
	}


	//FILE *fout = fopen_s("testVOL.vol", "wb");
	FILE *fout;
	fopen_s(&fout, volumeFile, "wb");
	static const char tag[] = "VOL";
	fwrite(tag, 1, 3, fout);
	static const unsigned char ver = 0x3;
	fwrite(&ver, sizeof(uint8_t), 1, fout);
	int data_format = 1;
	fwrite(&data_format, sizeof(int), 1, fout);

	// Write resolution
	fwrite(resol, sizeof(int), 3, fout);

	int ch = 1;
	fwrite(&ch, sizeof(int), 1, fout);

	// Write AABB
	//const float crop = 0.15; //needed for cropping stretching dataset
	//minAABB[0] = min_x - radius + crop*max_x, minAABB[1] = min_y - radius + crop*max_y, minAABB[2] = min_z - radius + crop*max_z; //use max as min is close to zero
	//maxAABB[0] = max_x + radius - crop*max_x, maxAABB[1] = max_y + radius - crop*max_y, maxAABB[2] = max_z + radius - crop*max_z;
	fwrite(minAABB, sizeof(float), 3, fout);
	fwrite(maxAABB, sizeof(float), 3, fout);

	// write voxel extent
	for (int i = 0; i < N; i++)
		fwrite(&(data[i]), sizeof(float), 1, fout);
	delete[] data;

	fclose(fout);

}


void step0_curveSetup(const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const int upsample) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** step0: up-sample the curve ***\n";
	const int num_of_cores = omp_get_num_procs();

	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;
		std::cout << "up-sample frame " << f << " started...\n";
		
		HermiteCurve curve_ds, curve_us;
		int seg_subdiv = 10;
		for (int y = yarn0; y < yarn1; ++y) {

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

			int seg_subdiv = 10;
			HermiteCurve curve_ds;
			curve_ds.init(curvefile_ds, normfile_ds, seg_subdiv);
			std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;

			// assign local frames for each point
			curve_ds.assign_twist(twistfile, all_pts, all_tang, all_norm, upsample);

			// write up-sampled centerline so later we can crop a segment out of it
			fout_cntr << vrtx << std::endl;
			fout_norm << vrtx << std::endl;
			for (int v = 0; v < vrtx; ++v) {
				fout_cntr << all_pts[v][0] << " " << all_pts[v][1] << " " << all_pts[v][2] << "\n";
				fout_norm << all_norm[v][0] << " " << all_norm[v][1] << " " << all_norm[v][2] << "\n";
			}
			fout_cntr.close();
			fout_norm.close();

			/* debug purpose */
			//std::ofstream fout_TNB("../data/TNB.txt");
			//for (int i = 0; i < vrtx_num; i++) {
			//	fout_TNB << all_pts[i][0] << " " << all_pts[i][1] << " " << all_pts[i][2] << " " <<
			//		all_tang[i][0] << " " << all_tang[i][1] << " " << all_tang[i][2] << " " <<
			//		all_norm[i][0] << " " << all_norm[i][1] << " " << all_norm[i][2] << "\n";
			//}

		}
	}
}

void step1_dg2local(const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** step1: Convert external-force to local coordinate ***\n";
	const int num_of_cores = omp_get_num_procs();

	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;

		std::cout << "Transfer DG to local for frame " << f << "\n";
		HermiteCurve curve;
		int seg_subdiv = 10;
		for (int y = yarn0; y < yarn1; ++y) {

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
			std::ifstream fin4(physical_world);
			assert(fin4.is_open() && "physical_world file wasn't found!\n");

			curve.init_norm(curvefile_us, normfile_us, seg_subdiv); /******* normals are already twisted *******/
			std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;
			curve.assign(all_pts, all_tang, all_norm);
			//curve.assign_twist(twistfile, all_pts, all_tang, all_norm, 1);
			assert(all_pts.size() == vrtx);
			//curve.init_principleNormal(curvefile, normfile, seg_subdiv);


			std::ifstream fin(physical_world);
			std::ofstream fout(physical_local);

			float S00, S01, S02, S10, S11, S12, S20, S21, S22;
			//float A0, A1, A2;
			//float B0, B1, B2;
			const int vrtx_num = vrtx;

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
				M << ez[0], ez[1], ez[2],
					ey[0], ey[1], ey[2],
					ex[0], ex[1], ex[2];

				//local = M*world*M.transpose();
				local = M*world;

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
}
void step1_shapematching(const char* yarnfile1, const char* configfile, const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const float stepSize) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** step1: Fitting phase ***\n";

	/* This yarn is what will be compressed (has flyaways) */
	Fiber::Yarn yarn;
	yarn.parse(configfile);
	yarn.setStepNum(vrtx);
	yarn.setStepSize(stepSize);
	yarn.simulate_ply_shuang();
	yarn.write_plys("test_ply.txt");
	const int K = yarn.getPlyNum();
	yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
	yarn.build("test_fly.txt", K);


	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;


		for (int y = yarn0; y < yarn1; ++y) {

			std::string tmp1 = "data/" + dataset + "/simul_frame_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* yarnfile2 = tmp1.c_str();
			std::string tmp3 = "input/" + dataset + "/matrix_S_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S = tmp3.c_str();
			std::string tmp6 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* normfile_us = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* curvefile_us = tmp7.c_str();
			std::string tmp5 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* twistfile = tmp5.c_str();
			std::string tmp10 = "input/" + dataset + "/physical_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local = tmp10.c_str();

			std::string tmp11 = "input/" + dataset + "/global_rot_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* global_rot = tmp11.c_str();

			std::ifstream fin1(yarnfile1);
			std::ifstream fin2(yarnfile2);

			assert(fin1.is_open() && "reference-yarn file wasn't found!\n");
			assert(fin2.is_open() && "compressed-yarn file wasn't found!\n");

			const int vrtx_num = yarn.getStepNum();

			const int upsample = 2;
			std::cout << "#### shapematching frame" << f << " yarn" << "  is started... \n";
			extractCompress_seg(configfile, yarnfile1, yarnfile2, "noNeed.txt", compress_S,
				curvefile_us, normfile_us, global_rot, yarn.getPlyNum(), vrtx_num);

			/*************************************************/
			std::string tmp8, tmp9;
#ifndef IMPROVED_FLYAWAYS
			tmp8 = "output/" + dataset + "/genYarn_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			tmp9 = "output/" + dataset + "/genYarn_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#else
			tmp8 = "output/" + dataset + "/genYarn_fly_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			tmp9 = "output/" + dataset + "/genYarn_wo_fly_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#endif

#if 1		

			//const char* outfile = tmp8.c_str();
			////// Procedural step
			//Fiber::Yarn yarn_compressed; //renew the yarn
			//yarn_compressed = yarn;
			//yarn_compressed.compress_yarn_A(compress_S, global_rot);
			//yarn_compressed.curve_yarn(curvefile_us, normfile_us);
			//yarn_compressed.write_yarn(outfile);

			/////*************************************************/

			const char* outfile_wo = tmp9.c_str();
			Fiber::Yarn yarn_wo; //renew the yarn
			yarn_wo = yarn;
			yarn_wo.rotate_yarn(global_rot);
			yarn_wo.curve_yarn(curvefile_us, normfile_us);
			yarn_wo.write_yarn(outfile_wo);
#endif

		}
	}
}

void step2_buildTrainData(const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const int isTrain, const int ws_ds, const float trimPercent, const int upsample) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** step 2: Build training data  *** \n";
	std::cout << " @@@@@@@@@@ " << dataset << " @@@@@@@@@@ \n";
	const int num_of_cores = omp_get_num_procs();

	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;

		int seg_subdiv = 10;
		std::cout << "Build training data for frame " << f << " started ... \n";

		for (int y = yarn0; y < yarn1; ++y) {

			std::string tmp3 = "input/" + dataset + "/physicalParam/physical_" + std::to_string(f) + "_" + std::to_string(y) + "_world.txt";
			const char* physical_world = tmp3.c_str();
			std::string tmp4 = "input/" + dataset + "/NN/trainX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg = tmp4.c_str();
			std::string tmp5 = "input/" + dataset + "/matrix_S_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S = tmp5.c_str();
			std::string tmp6 = "input/" + dataset + "/NN/trainY_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S_seg = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/NN/angles_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* angles = tmp7.c_str();

			std::string tmp8 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* twistfile = tmp8.c_str();
			std::string tmp9 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* curvefile_us = tmp9.c_str();
			std::string tmp10 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* normfile_us = tmp10.c_str();

			std::string tmp11 = "input/" + dataset + "/NN/testX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_test = tmp11.c_str();

			std::ifstream fin3(physical_world);
			assert(fin3.is_open() && "physical_world file wasn't found!\n");

			const int vrtx_num = vrtx;

			std::ofstream fout_testX(physical_local_seg_test);
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
			const int window_size = ws_ds * upsample;
			


			const int window_num = ((vrtx_num - window_size + 1) - 2 * ignorPlanes);

			std::vector<int> indices;
			for (int w = ignorPlanes; w < (vrtx_num - window_size + 1) - ignorPlanes; w = w + upsample)
				indices.push_back(w);
			std::vector <Eigen::VectorXd> store_NN_test_total(indices.size());
			std::vector <float> angle_total(indices.size());

			const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
			for (int omp_i = 0; omp_i < static_cast<int>(indices.size()); ++omp_i) {
				const int w = indices[omp_i];

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

				std::vector<Eigen::Matrix3f> all_dg_seg;
				for (int d = start; d <= end; d++) all_dg_seg.push_back(all_dg[d]);

				std::vector<Eigen::Matrix3f> all_local_dg_seg, all_local_dg_seg_flip_norm, all_local_dg_seg_flip_tang, all_local_dg_seg_flip_both;
				transfer_dg_2local(all_tang_seg, all_norm_seg, all_dg_seg, all_local_dg_seg, 0);
				if (isTrain) {
					transfer_dg_2local(all_tang_seg, all_norm_seg, all_dg_seg, all_local_dg_seg_flip_norm, 1); //augment the data by including the flipped normals
					transfer_dg_2local(all_tang_seg, all_norm_seg, all_dg_seg, all_local_dg_seg_flip_tang, 2); //augment the data by including the flipped tangents
					transfer_dg_2local(all_tang_seg, all_norm_seg, all_dg_seg, all_local_dg_seg_flip_both, 3); //augment the data by including the flipped tangents
				}


				if (isTrain) {
					for (int d = 0; d < window_size; (d += upsample)) {
						Eigen::Matrix3f local_dg = all_local_dg_seg[d];
						fout_trainX << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
							local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
							local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
					}
					fout_trainX << "\n";
				}

				const int v_yarn = ceil((start + end) / 2.0);
				const int v_seg = ceil((end - start) / 2.0);
				assert(v_yarn == v_seg + start && "index for full yarn must be equal to index segment added with starting vertex");

				Eigen::Vector3d norm1 = all_norm[v_yarn];
				Eigen::Vector3d norm2 = all_norm_seg[v_seg];
				Eigen::Vector3d tang = all_tang[v_yarn];
				assert((all_tang[v_yarn] - all_tang_seg[v_seg]).norm() < eps && "tangents for both frames must be similar!\n");
				float angle = get_angle(norm1, norm2, tang);
				/****************************/
				angle_total[omp_i] = angle;
				//fout_angle << angle << "\n";

				if (isTrain) {
					/************* write train y **************/
					if (isTrain) {
						Eigen::Matrix2f S_local;
						rotate_S_2local(all_S[v_yarn], S_local, angle, 0);
						fout_trainY << S_local(0, 0) << " " << S_local(0, 1) << " " << S_local(1, 0) << " " << S_local(1, 1) << "\n";
					}
				}

				// d+=upsample because we skip the upsampled DG's and we only care about the simulated ones
				
				/******** write test data *******/
				///Eigen::VectorXd store_NN_test(int(window_size * 9)); //9 is size of DG
				Eigen::VectorXd store_NN_test(int(ws_ds * 9)); //9 is size of DG
				int p = 0;
				for (int d = 0; d < window_size; (d += upsample)) {

					Eigen::Matrix3f local_dg = all_local_dg_seg[d];

					store_NN_test(p * 9 + 0) = local_dg(0, 0);
					store_NN_test(p * 9 + 1) = local_dg(0, 1);
					store_NN_test(p * 9 + 2) = local_dg(0, 2);
					store_NN_test(p * 9 + 3) = local_dg(1, 0);
					store_NN_test(p * 9 + 4) = local_dg(1, 1);
					store_NN_test(p * 9 + 5) = local_dg(1, 2);
					store_NN_test(p * 9 + 6) = local_dg(2, 0);
					store_NN_test(p * 9 + 7) = local_dg(2, 1);
					store_NN_test(p * 9 + 8) = local_dg(2, 2);
					p++;
				}
				store_NN_test_total[omp_i] = store_NN_test;


				///******** write test data *******/
				//for (int d = 0; d < window_size; (d += upsample)) {

				//	Eigen::Matrix3f local_dg = all_local_dg_seg[d];
				//	fout_testX << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
				//		local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
				//		local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
				//}
				//fout_testX << "\n";

				if (isTrain) {
#if 1 /* augment the training */
					
					//****** augment by flipping normals ******
					for (int d = 0; d < window_size; (d += upsample)) {

						Eigen::Matrix3f local_dg = all_local_dg_seg_flip_norm[d];
						fout_trainX << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
							local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
							local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
					}
					fout_trainX << "\n";

					if (isTrain) {
						Eigen::Matrix2f S_local;
						rotate_S_2local(all_S[v_yarn], S_local, angle, 1);
						fout_trainY << S_local(0, 0) << " " << S_local(0, 1) << " " << S_local(1, 0) << " " << S_local(1, 1) << "\n";
					}
					//****** augment by flipping tangents ******
					for (int d = 0; d < window_size; (d += upsample)) {
						Eigen::Matrix3f local_dg = all_local_dg_seg_flip_tang[d];
						fout_trainX << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
							local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
							local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
					}
					fout_trainX << "\n";


					if (isTrain) {
						Eigen::Matrix2f S_local;
						rotate_S_2local(all_S[v_yarn], S_local, angle, 2);
						fout_trainY << S_local(0, 0) << " " << S_local(0, 1) << " " << S_local(1, 0) << " " << S_local(1, 1) << "\n";
					}
					//****** augment by flipping tangents and normals ******
					for (int d = 0; d < window_size; (d += upsample)) {
						Eigen::Matrix3f local_dg = all_local_dg_seg_flip_both[d];
						fout_trainX << local_dg(0, 0) << " " << local_dg(0, 1) << " " << local_dg(0, 2) << " " <<
							local_dg(1, 0) << " " << local_dg(1, 1) << " " << local_dg(1, 2) << " " <<
							local_dg(2, 0) << " " << local_dg(2, 1) << " " << local_dg(2, 2) << " ";
					}
					fout_trainX << "\n";


					if (isTrain) {
						Eigen::Matrix2f S_local;
						rotate_S_2local(all_S[v_yarn], S_local, angle, 3);
						fout_trainY << S_local(0, 0) << " " << S_local(0, 1) << " " << S_local(1, 0) << " " << S_local(1, 1) << "\n";
					}
#endif
				}
			}

			for (int l = 0; l < store_NN_test_total.size(); l++) {
				//for (int p = 0; p < int(window_size * 9); p++) { 
				// this was writing all window-size although we want the downsampled Dg's which are the actual 
				// simulated ones. Have to sanity check if I mess up other parts
				for (int p = 0; p < int(ws_ds * 9); p++) {
					fout_testX << store_NN_test_total[l](p) << " ";
				}
				fout_testX << "\n";
			}
			for (int l = 0; l < store_NN_test_total.size(); l++) {
				fout_angle << angle_total[l] << "\n";
			}

			fout_testX.close();
			fout_trainX.close();
			fout_trainY.close();
			fout_angle.close();
		}

	}
}


void step3_temporalTrainData(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const int isTrain) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** step 2: Add temporal support to training data  *** \n";
	std::cout << " @@@@@@@@@@ " << dataset << " @@@@@@@@@@ \n";

	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;
		for (int y = yarn0; y < yarn1; ++y) {

			
			// Padding for the first frame 
			int preFrame = f - skipFactor;
			int postFrame = f + skipFactor;
			if (i == 0)
				preFrame = f;
			else if (i == total_frame - 1)
				postFrame = f;

			std::string tmp0 = "input/" + dataset + "/NN/trainX_" + std::to_string(preFrame) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_pre = tmp0.c_str();
			std::string tmp1 = "input/" + dataset + "/NN/testX_" + std::to_string(preFrame) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_test_pre = tmp1.c_str();

			std::string tmp2 = "input/" + dataset + "/NN/trainX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg = tmp2.c_str();
			std::string tmp3 = "input/" + dataset + "/NN/testX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_test = tmp3.c_str();

			std::string tmp4 = "input/" + dataset + "/NN/trainX_" + std::to_string(postFrame) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_post = tmp4.c_str();
			std::string tmp5 = "input/" + dataset + "/NN/testX_" + std::to_string(postFrame) + "_" + std::to_string(y) + ".txt";
			const char* physical_local_seg_test_post = tmp5.c_str();

			std::string tmp6 = "input/" + dataset + "/NN/trainX_" + std::to_string(f) + "_" + std::to_string(y) + "_temporal.txt";
			const char* physical_local_seg_temporal = tmp6.c_str();
			std::string tmp7 = "input/" + dataset + "/NN/testX_" + std::to_string(f) + "_" + std::to_string(y) + "_temporal.txt";
			const char* physical_local_seg_test_temporal = tmp7.c_str();


			if (isTrain) {
				std::ifstream finTrain_pre(physical_local_seg_pre);
				assert(finTrain_pre.is_open() && "trainx pre wasn't found!\n");
				std::ifstream finTrain(physical_local_seg);
				assert(finTrain.is_open() && "trainx wasn't found!\n");
				std::ifstream finTrain_post(physical_local_seg_post);
				assert(finTrain_post.is_open() && "trainx post wasn't found!\n");

				std::ofstream foutTrain(physical_local_seg_temporal);


				for (int i = 0; finTrain.eof() != true; i++) {
					std::string content_trainX = "";
					std::string content_trainX_pre = "";
					std::string content_trainX_post = "";

					std::getline(finTrain_pre, content_trainX_pre);
					std::getline(finTrain, content_trainX);
					std::getline(finTrain_post, content_trainX_post);

					foutTrain << content_trainX_pre << " " << content_trainX << " " << content_trainX_post << std::endl;

				}
				//content_trainX_temporal.erase(content_trainX_temporal.end() - 1);     // erase last character
				finTrain_pre.close();
				finTrain.close();
				finTrain_post.close();

				// write the updated training data

				foutTrain.close();
			}
			std::ifstream finTest_pre(physical_local_seg_test_pre);
			assert(finTest_pre.is_open() && "testx pre wasn't found!\n");
			std::ifstream finTest(physical_local_seg_test);
			assert(finTest.is_open() && "testx wasn't found!\n");
			std::ifstream finTest_post(physical_local_seg_test_post);
			assert(finTest_post.is_open() && "testx post wasn't found!\n");

			std::ofstream foutTest(physical_local_seg_test_temporal);


			for (int i = 0; finTest.eof() != true; i++) {
				std::string content_testX = "";
				std::string content_testX_pre = "";
				std::string content_testX_post = "";

				std::getline(finTest_pre, content_testX_pre);
				std::getline(finTest, content_testX);
				std::getline(finTest_post, content_testX_post);

				foutTest << content_testX_pre << " " << content_testX << " " << content_testX_post << std::endl ;

			}
			//content_trainX_temporal.erase(content_trainX_temporal.end() - 1);     // erase last character
			finTest_pre.close();
			finTest.close();
			finTest_post.close();

			// write the updated training data

			foutTest.close();
		}
	}
}


void step3_appendTraining(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** step 3: Append training-data for all frames *** \n";

	std::string tmp0 = "input/" + dataset + "/NN/trainX_all.txt";
	const char* all_trainX = tmp0.c_str();
	std::string tmp1 = "input/" + dataset + "/NN/trainY_all.txt";
	const char* all_trainY = tmp1.c_str();
	std::ofstream fout_trainX_all(all_trainX);
	std::ofstream fout_trainY_all(all_trainY);
	std::string content_trainX = "";
	std::string content_trainY = "";

	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;

		int seg_subdiv = 10;
		for (int y = yarn0; y < yarn1; ++y) {
			std::string tmp4 = "input/" + dataset + "/NN/trainX_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* trainX = tmp4.c_str();
			std::string tmp6 = "input/" + dataset + "/NN/trainY_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* trainY = tmp6.c_str();


			std::ifstream finX(trainX);
			assert(finX.is_open() && "trainX file wasn't found!\n");
			std::ifstream finY(trainY);
			assert(finY.is_open() && "trainY file wasn't found!\n");

			int i;
			for (i = 0; finX.eof() != true; i++) // get content of infile
				content_trainX += finX.get();
			i--;
			content_trainX.erase(content_trainX.end() - 1);     // erase last character
			finX.close();


			// trainY
			int j;
			for (j = 0; finY.eof() != true; j++) // get content of infile
				content_trainY += finY.get();
			j--;
			content_trainY.erase(content_trainY.end() - 1);     // erase last character
			finY.close();


		}
	}
	fout_trainX_all << content_trainX;                 // output
	fout_trainY_all << content_trainY;                 // output
	fout_trainX_all.close();
	fout_trainY_all.close();

}

void step4_NN_output(const char* configfile, const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const int isTrain, const int isCompress, const float stepSize) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** Testing-NN phase ***\n";
	std::cout << " @@@@@@@@@@ " << dataset << " @@@@@@@@@@ \n";
	const int num_of_cores = omp_get_num_procs();

	//// Procedural step
	/* This yarn is what will be compressed (has flyaways) */
	Fiber::Yarn yarn;
	yarn.parse(configfile);
	yarn.setStepNum(vrtx);
	yarn.setStepSize(stepSize);
	yarn.simulate_ply_shuang();
	yarn.write_plys("test_ply.txt");
	const int K = yarn.getPlyNum();
	yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
	yarn.build("test_fly.txt", K);

	const int total_frame = (frame1 - frame0) / skipFactor + 1; 
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;

		std::cout << "Generate fibers for frame " << f << " started... \n";

		for (int y = yarn0; y < yarn1; ++y) {

			std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* curvefile_us = tmp1.c_str();
			std::string tmp2 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
			const char* normfile_us = tmp2.c_str();
			std::string tmp6 = "input/" + dataset + "/NN/testY_NN_full_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* compress_S = tmp6.c_str();

			std::string tmp11 = "input/" + dataset + "/global_rot_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
			const char* global_rot = tmp11.c_str();

			std::cout << compress_S << std::endl;
			std::ifstream fin2(compress_S);
			if (isTrain)
				assert(fin2.is_open() && "compress_S_NN file wasn't found!\n");
			std::ifstream fin3(curvefile_us);
			assert(fin3.is_open() && "curvefile file wasn't found!\n");
			std::ifstream fin4(normfile_us);
			assert(fin4.is_open() && "normfile file wasn't found!\n");

			///*******  write the yarn ******/
			std::string tmp3;
			if (isCompress) {
#ifndef IMPROVED_FLYAWAYS
				tmp3 = "output/" + dataset + "/genYarn_NN_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#else
				tmp3 = "output/" + dataset + "/genYarn_NN_fly_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#endif
			}
			else {
#ifndef IMPROVED_FLYAWAYS
				tmp3 = "output/" + dataset + "/genYarn_NN_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#else
				tmp3 = "output/" + dataset + "/genYarn_NN_fly_wo_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
#endif
			}

			const char* outfile = tmp3.c_str();
			////// Procedural step
			Fiber::Yarn yarn_compressed; //renew the yarn
			yarn_compressed = yarn;
			if (isTrain)
				yarn_compressed.compress_yarn_A(compress_S, global_rot);
			else {
				if (isCompress)
					yarn_compressed.compress_yarn_A(compress_S);
			}
			yarn_compressed.curve_yarn(curvefile_us, normfile_us);
			yarn_compressed.write_yarn(outfile);

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

void step5_createVOL(int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset, const int resol_x, const int resol_y, const int resol_z, const float radius) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** Create volume phase ***\n";
	std::cout << " @@@@@@@@@@ " << dataset << " @@@@@@@@@@ \n";

	const int num_of_cores = omp_get_num_procs();

	const int total_frame = (frame1 - frame0) / skipFactor + 1;
	for (int i = 0; i < total_frame; i++) {
		int f = frame0 + i * skipFactor;

		std::string tmp = "output/" + dataset + "/volume_" + std::to_string(f) + ".vol";
		const char* volfile_us = tmp.c_str();
		std::cout << frame0 << " " << f << " " << volfile_us << " generation is started... \n";
		int curr_frame = f;
		writeVol(dataset, curr_frame, yarn0, yarn1, resol_x, resol_y, resol_z, radius, volfile_us);
	}
}

void full_pipeline(const char* yarnfile1, const char* configfile, const int vrtx, int skipFactor, int frame0, int frame1, int yarn0, int yarn1, std::string &dataset,
	const int isTrain, const int window_size, const float trimPercent, const int upsample, const float stepSize) {
	std::cout << " @@@@@@@@@@ " << dataset << " @@@@@@@@@@ \n";
	step0_curveSetup(vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, upsample);
	step1_dg2local(vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
	if (isTrain)
		step1_shapematching(yarnfile1, configfile, vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, stepSize);
	step2_buildTrainData(vrtx, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent, upsample);
	//step3_temporalTrainData(skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain);
	if (isTrain)
		step3_appendTraining(skipFactor, frame0, frame1, yarn0, yarn1, dataset);
}


void upsampleMatrix(const char* filename, const int upsample, const char* filename_us, const int offset) {
	// file contains 4 values in each row
	std::ifstream fin(filename);
	if (!fin.is_open()) std::cout << filename << " not fount! \n";
	assert(fin.is_open());

	Eigen::Vector4f value;
	std::vector<Eigen::Vector4f> values;
	while (1) {
		fin >> value[0] >> value[1] >> value[2] >> value[3];
		if (fin.eof()) break;
		values.push_back(value);
	}
	const int N = values.size();

	// Let's interpolate for upsampled vector
	std::vector<Eigen::Vector4f> values_us;
	for (int j = 0; j < N; j++) {
		for (int s = 0; s < upsample; s++) {
			if (j == N - 1) {
				values_us.push_back(values[j]);
			}
			else {
				const float w = float(s) / float(upsample);
				Eigen::Vector4f value_us = (1.0 - w) * values[j] + w * values[j + 1];
				//std::cout << values[j] << " \n" << values[j + 1] << " \n" << value_us << std::endl << std::endl;
				values_us.push_back(value_us);
			}
		}
	}


	assert(values_us.size() == N*upsample);
	std::ofstream fout(filename_us);
	//fout << N*upsample << "\n";
	for (int i = 0; i < N*upsample; i++) {
		if (offset) { // rotation_90 = [0 -1 ; 1 0] so [a b; c d][0 -1; 1 0] = [b -a; d -c]
			fout << values_us[i][1] << " " << -1.0*values_us[i][0] << " " << values_us[i][3] << " " << -1.0*values_us[i][2] << "\n";
			continue;
		}
		fout << values_us[i][0] << " " << values_us[i][1] << " " << values_us[i][2] << " " << values_us[i][3] << "\n";
	}
	fout.close();
}

void upsampleValue(const char* filename, const int upsample, const char* filename_us) {
	// file starts with #N and contains 1 value in a row 
	std::ifstream fin(filename);
	if (!fin.is_open()) std::cout << filename;
	assert(fin.is_open());

	int N;
	float value;
	std::vector<float> values;
	fin >> N;
	for (int i = 0; i < N; i++) {
		fin >> value;
		values.push_back(value);
	}

	// Let's interpolate for upsampled vector
	std::vector<float> values_us;
	for (int j = 0; j < N; j++) {
		for (int s = 0; s < upsample; s++) {
			if (j == N - 1) {
				values_us.push_back(values[j]);
			}
			else {
				const float w = float(s) / float(upsample);
				float value_us = (1.0 - w) * values[j] + w * values[j + 1];
				//std::cout << values[j] << " " << values[j+1] << " " << value_us << std::endl;

				values_us.push_back(value_us);
			}
		}
	}

	assert(values_us.size() == N*upsample);
	std::ofstream fout(filename_us);
	fout << N*upsample << "\n";
	for (int i = 0; i < N*upsample; i++) {
		fout << values_us[i] << "\n";
	}
	fout.close();

}

void upsampleCurve(const int vrtx, const int upsample, const char* curvefile_ds, const char* curvefile_us, const char* normfile_ds, const char* normfile_us, const char* twistfile) {
	std::cout << "\n**************************************************\n";
	std::cout << "*** step0: up-sample the curve for stretched yarn ***\n";

	HermiteCurve curve_ds, curve_us;
	int seg_subdiv = 10;


	std::ofstream fout_cntr(curvefile_us);
	std::ofstream fout_norm(normfile_us);

	curve_ds.init_norm(curvefile_ds, normfile_ds, seg_subdiv);
	std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;

	// assign local frames for each point
	curve_ds.assign_twist(twistfile, all_pts, all_tang, all_norm, upsample);

	// write up-sampled centerline so later we can crop a segment out of it
	fout_cntr << vrtx << std::endl;
	fout_norm << vrtx << std::endl;
	for (int v = 0; v < vrtx; ++v) {
		fout_cntr << all_pts[v][0] << " " << all_pts[v][1] << " " << all_pts[v][2] << "\n";
		fout_norm << all_norm[v][0] << " " << all_norm[v][1] << " " << all_norm[v][2] << "\n";
	}
	fout_cntr.close();
	fout_norm.close();
}

void upsample_stretched(const char* configfile, const int vrtx, const int f, std::string &dataset, const int upsample) {

	int y = 0;
	HermiteCurve curve_ds, curve_us;
	int seg_subdiv = 10;

	std::string tmp1 = "input/" + dataset + "/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
	const char* curvefile_ds = tmp1.c_str();
	std::string tmp2 = "input/" + dataset + "/stretch/centerYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us_us.txt";
	const char* curvefile_us = tmp2.c_str();
	std::string tmp3 = "input/" + dataset + "/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
	const char* normfile_ds = tmp3.c_str();
	std::string tmp4 = "input/" + dataset + "/stretch/normYarn_" + std::to_string(f) + "_" + std::to_string(y) + "_us_us.txt";
	const char* normfile_us = tmp4.c_str();
	std::string tmp5 = "input/" + dataset + "/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
	const char* twistfile_ds = tmp5.c_str();
	std::string tmp6 = "input/" + dataset + "/stretch/twist_" + std::to_string(f) + "_" + std::to_string(y) + "_us_us.txt";
	const char* twistfile = tmp6.c_str();
	std::string tmp7 = "input/" + dataset + "/NN/testY_NN_full_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
	const char* compress_S_ds = tmp7.c_str();
	std::string tmp8 = "input/" + dataset + "/stretch/testY_NN_full_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
	const char* compress_S = tmp8.c_str();
	std::string tmp9 = "input/" + dataset + "/global_rot_" + std::to_string(f) + "_" + std::to_string(y) + ".txt";
	const char* global_rot_ds = tmp9.c_str();
	std::string tmp10 = "input/" + dataset + "/stretch/global_rot_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
	const char* global_rot = tmp10.c_str();
	std::string tmp11 = "output/" + dataset + "/genYarn_NN_" + std::to_string(f) + "_" + std::to_string(y) + "_us.txt";
	const char* outfile = tmp11.c_str();

	upsampleValue(twistfile_ds, upsample, twistfile);
	upsampleMatrix(compress_S_ds, upsample, compress_S, 0);
	upsampleMatrix(global_rot_ds, upsample, global_rot, 0);
	upsampleCurve(vrtx, upsample, curvefile_ds, curvefile_us, normfile_ds, normfile_us, twistfile);


	//// Procedural step
	Fiber::Yarn yarn;
	yarn.parse(configfile);
	yarn.setStepNum(vrtx);
	const float stepSize = yarn.getStepSize() / float(upsample);
	yarn.setStepSize(stepSize);
	yarn.simulate_ply_shuang();
	yarn.write_plys("test_ply.txt");
	const int K = yarn.getPlyNum();
	yarn.roll_plys(K, "test_ply.txt", "test_fly.txt");
	yarn.build("test_fly.txt", K);

	Fiber::Yarn yarn_compressed; //renew the yarn
	yarn_compressed = yarn;
	yarn_compressed.compress_yarn_A(compress_S, global_rot);
	yarn_compressed.curve_yarn(curvefile_us, normfile_us);
	yarn_compressed.write_yarn(outfile);

}

#if 0
void buildTraining(const char* curvefile_ds, const char* normfile_ds, const char* physical_world, const char* compress_S, Fiber::Yarn &yarn, const float trimPercent,
	const int window_size, const char* curvefile_us, const char* angles, const char* physical_local_window, const char* physical_local_window_test, const char* compress_S_window,
	std::ofstream &fout_trainX_all, std::ofstream &fout_trainY_all, int isTrain) {

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

# if 1
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
void buildTraning_all(Fiber::Yarn &yarn, int skipFactor, int frame0, int frame1, int yarnNum, std::string &dataset, const int window_size, const float trimPercent, const int isTrain) {

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
#endif

#if 0
void constFitting_compParam(const std::vector<Ellipse> &ellipses, const std::vector<float> &theta_R,
	const float trimPercent, Fiber::Yarn::Compress &compress) {

	const int ignorPlanes = trimPercent * ellipses.size();
	Point2DVector points_R1;
	Point2DVector points_R2;
	Point2DVector points_theta;
	Point2DVector points_rot;

	for (int i = ignorPlanes; i < ellipses.size() - ignorPlanes; ++i) {
		Eigen::Vector2d point;
		point(0) = i - ignorPlanes;
		point(1) = ellipses[i].longR;
		points_R1.push_back(point);

		point(1) = ellipses[i].shortR;
		points_R2.push_back(point);

		point(1) = ellipses[i].angle;
		points_theta.push_back(point);

		point(1) = theta_R[i];
		points_rot.push_back(point);
	}
	Eigen::VectorXd x(1);
	x.fill(0.f);
	MyFunctorNumericalDiff functor;
	functor.Points = points_R1;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);
	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(x);
	compress.ellipse_long = x(0);

	functor.Points = points_R2;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm2(functor);
	status = lm2.minimize(x);
	compress.ellipse_short = x(0);

	functor.Points = points_theta;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm3(functor);
	status = lm3.minimize(x);
	compress.ellipse_theta = x(0);

	functor.Points = points_rot;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm4(functor);
	status = lm4.minimize(x);
	compress.rotation = x(0);

	//std::cout << "fitted compression params: " << compress.ellipse_long << " " << compress.ellipse_short << " " << compress.ellipse_theta << " " << compress.rotation << "\n";
}

void sinFitting_curve(const char* curveFile, const float trimPercent, Fiber::Yarn::CenterLine &curve) {
	//read the curve points
	std::vector<float> curvePnts;
	std::ifstream fin;
	if (curveFile != NULL)
		fin.open(curveFile);
	std::string line;
	std::getline(fin, line);
	const int plane_num = atof(line.c_str());
	for (int i = 0; i < plane_num; ++i) {
		std::getline(fin, line);
		float pnt = atof(line.c_str());
		curvePnts.push_back(pnt);
	}

	const int ignorPlanes = trimPercent * curvePnts.size();
	Point2DVector points;

	for (int i = ignorPlanes; i < curvePnts.size() - ignorPlanes; ++i) {
		Eigen::Vector2d point;
		point(0) = i - ignorPlanes;
		point(1) = curvePnts[i];
		points.push_back(point);
	}
	const int period = curvePnts.size() - 2 * ignorPlanes;
	Eigen::VectorXd x(3);
	x(0) = 1.f;
	x(1) = 2.f * pi / static_cast<float>(period);
	x(2) = 0.f;


	MyFunctorNumericalDiff functor;
	functor.Points = points;
	Eigen::LevenbergMarquardt<MyFunctorNumericalDiff> lm(functor);
	Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(x);
	curve.a = x(0);
	//curve.b = x(1);
	curve.b = 2.f * pi / static_cast<float>(period); //TODO: pass this during fitting
	curve.c = x(2);

	//std::cout << "curve params: " << curve.a << " " << curve.b << " compared to: " << 2.f * pi / static_cast<float>(period) << " " << curve.c << std::endl;

	//for debug: 
	std::ofstream fout("sineCurve.txt");
	fout << points.size() << std::endl;
	for (int i = 0; i < points.size(); ++i) {
		float y = curve.a * cos(curve.b*i) + curve.c;
		fout << y << std::endl;
	}
	fout.close();
}
#endif 

#if 0
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
#endif 

#if 0
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
#endif

#if 0
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
#endif


#if 0
void interpolate(const std::vector<float> vals, const int interval, std::vector<float> newVals) {

}

void appendCompress_yarn(const std::vector<Fiber::Yarn::Compress> &compress_segs, const int seg_vrtx,
	const int yarn_vrtx, const char* compressFile) {
	const int seg_num = compress_segs.size();
	assert(seg_num * seg_vrtx < yarn_vrtx);
	std::vector<Fiber::Yarn::Compress> allCompress;
	allCompress.resize(yarn_vrtx);

	const int seg_vrtx_new = yarn_vrtx / seg_num;
	//interpolate old segment to new one
	int indx = 0;
	for (int i = 0; i < seg_num; ++i) {
		for (int v = 0; v < seg_vrtx_new; ++v) {
			indx = i*seg_vrtx_new + v;
			allCompress[indx].ellipse_long = compress_segs[i].ellipse_long;
			allCompress[indx].ellipse_short = compress_segs[i].ellipse_short;
			allCompress[indx].ellipse_theta = compress_segs[i].ellipse_theta;
			allCompress[indx].rotation = compress_segs[i].rotation;
		}
	}
	//mirror the values for the leftover vertices  (if the end of yarn hasn't been assign)
	int c = 1;
	for (int l = indx + 1; l < yarn_vrtx; ++l) {
		allCompress[l].ellipse_long = allCompress[indx - c].ellipse_long;
		allCompress[l].ellipse_short = allCompress[indx - c].ellipse_short;
		allCompress[l].ellipse_theta = allCompress[indx - c].ellipse_theta;
		allCompress[l].rotation = allCompress[indx - c].rotation;
		c++;
	}

	//write to file
	FILE *fout;
	//write ellipses to file for testing
	if (fopen_s(&fout, compressFile, "wt") == 0) {
		fprintf_s(fout, "%d \n", yarn_vrtx);
		for (int i = 0; i < yarn_vrtx; ++i) {
			fprintf_s(fout, "%.4f %.4f %.4f %.4f\n", allCompress[i].ellipse_long, allCompress[i].ellipse_short, allCompress[i].ellipse_theta, allCompress[i].rotation);
		}
		fclose(fout);
	}
}

void appendCenter_yarn(const std::vector<Fiber::Yarn::CenterLine> &centerlines, const int seg_vrtx, const int yarn_vrtx, const char* curveFile) {
	const int seg_num = centerlines.size();
	assert(seg_num * seg_vrtx < yarn_vrtx);
	std::vector<Fiber::Yarn::CenterLine> allCenterlines;
	allCenterlines.resize(yarn_vrtx);
	const int seg_vrtx_new = yarn_vrtx / seg_num;
	//interpolate old segment to new one
	int indx = 0;
	for (int i = 0; i < seg_num; ++i) {
		for (int v = 0; v < seg_vrtx_new; ++v) {
			indx = i*seg_vrtx_new + v;
			allCenterlines[indx].a = centerlines[i].a;
			allCenterlines[indx].b = centerlines[i].b;
			allCenterlines[indx].c = centerlines[i].c;
		}
	}
	//mirror the values for the leftover vertices  (if the end of yarn hasn't been assign)
	int c = 1;
	for (int l = indx + 1; l < yarn_vrtx; ++l) {
		allCenterlines[l].a = centerlines[indx - c].a;
		allCenterlines[l].b = centerlines[indx - c].b;
		allCenterlines[l].c = centerlines[indx - c].c;
		c++;
	}

	//write to file
	FILE *fout;
	//write ellipses to file for testing
	if (fopen_s(&fout, curveFile, "wt") == 0) {
		fprintf_s(fout, "%d \n", yarn_vrtx);
		for (int i = 0; i < yarn_vrtx; ++i) {
			int v = i%seg_vrtx_new;
			float b = 2 * pi / (seg_vrtx_new); //stretch the curve to fit the new length
			float y = allCenterlines[i].a* sin(b * v) + allCenterlines[i].c;
			fprintf_s(fout, "0.0 %.4f \n", y); //TODO: this should be in 3D
											   //fprintf_s(fout, "%.4f %.4f %.4f\n", allCenterlines[i].a, allCenterlines[i].b, allCenterlines[i].c);
		}
		fclose(fout);
	}
}
#endif
