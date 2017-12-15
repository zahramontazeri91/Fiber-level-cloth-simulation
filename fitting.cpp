#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

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
	for (int l = indx+1; l < yarn_vrtx; ++l) {
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
			allCenterlines[indx].d = centerlines[i].d;
		}
	}
	//mirror the values for the leftover vertices  (if the end of yarn hasn't been assign)
	int c = 1;
	for (int l = indx + 1; l < yarn_vrtx; ++l) {
		allCenterlines[l].a = centerlines[indx - c].a;
		allCenterlines[l].b = centerlines[indx - c].b;
		allCenterlines[l].c = centerlines[indx - c].c;
		allCenterlines[l].d = centerlines[indx - c].d;
		c++;
	}

	//write to file
	FILE *fout;
	//write ellipses to file for testing
	if (fopen_s(&fout, curveFile, "wt") == 0) {
		fprintf_s(fout, "%d \n", yarn_vrtx);
		for (int i = 0; i < yarn_vrtx; ++i) {
			fprintf_s(fout, "%.4f %.4f %.4f %.4f\n", allCenterlines[i].a, allCenterlines[i].b, allCenterlines[i].c, allCenterlines[i].d);
		}
		fclose(fout);
	}
}

void extractCompress_seg(const char* yarnfile1, const char* yarnfile2, const char* compressFile, 
	const int ply_num, const int vrtx_num)
{
	const int n = vrtx_num;

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
	//std::cout << "\n";
	//std::cout << "Step 1: Shape matching between reference and compressed yarn.\n";
	std::vector<Ellipse> ellipses;
	std::vector<float> all_theta_R;
	cs.yarnShapeMatches(pnts_trans, pnts_ref, ellipses, all_theta_R);

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

	
	FILE *fout;
	if (fopen_s(&fout, compressFile, "wt") == 0) {
		fprintf_s(fout, "%d \n", ellipses.size());
		for (int i = 0; i < ellipses.size(); ++i) {
			//fprintf_s(fout, "%.6f %.6f %.6f ", validEllipses[i].longR, validEllipses[i].shortR, validEllipses[i].angle);
			//fprintf_s(fout, "%.6f %.6f %.6f \n", all_theta_R[i], all_T[i](0, 0), all_T[i](1, 0));
			fprintf_s(fout, "%.6f %.6f %.6f %.6f \n", ellipses[i].longR, ellipses[i].shortR, ellipses[i].angle, all_theta_R[i]);
		}
		fclose(fout);
	}
	std::cout << "Compression parameters are successfully written to the file!\n";

	//for debug: visualization
	const char* L2File = "../data/L2.txt";
	const char* refFile = "../data/allCrossSection2D_ref.txt";
	const char* deformedRefFile = "../data/allCrossSection2D_deformedRef.txt";
	const char* deformedFile = "../data/allCrossSection2D_deformed.txt";

	plotIntersections(pnts_ref, refFile);
	std::vector<yarnIntersect2D> ref_deformed;
	deformRef(pnts_ref, ref_deformed, ellipses, all_theta_R);
	plotIntersections(ref_deformed, deformedRefFile);
	plotIntersections(pnts_trans, deformedFile);

	std::vector<float> L2;
	L2norm(ref_deformed, pnts_trans, L2, L2File);
}

void regularize_seg() {
	//optimizeEllipses(const std::vector<Ellipse> &ellipses, const std::vector<float> &theta_R,
		//std::vector<Ellipse> &new_ellipses, std::vector<float> &new_theta_R);
}
void L2norm(const std::vector<yarnIntersect2D> &its_deform, const std::vector<yarnIntersect2D> &its_trans, std::vector<float> &L2, const char* filename) {
	assert(its_deform.size() == its_trans.size());
	FILE *fout;

	if (fopen_s(&fout, filename, "wt") == 0) {
		fprintf_s(fout, "%d \n", its_deform.size());
		for (int i = 0; i < its_deform.size(); ++i) { //number of planes
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

void plotIntersections(const std::vector<yarnIntersect2D> &its, const char* filename) {
	FILE *fout;
	// write the plycenters
	if (fopen_s(&fout, filename, "wt") == 0) {
		const int ignorPlanes = 0.2 * its.size(); // crop the first and last 10% of the yarn

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