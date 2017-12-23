#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"
#include "curveFitting.h"

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

void extractCompress_seg(const char* yarnfile1, const char* yarnfile2, const char* compressFile_vrtx, const char* compressFile_seg,
	const char* curveFile, const int ply_num, const int vrtx_num, std::vector<Ellipse> &new_ellipses, std::vector<float> &new_theta_R)
{
	const int n = vrtx_num;

	
	Fiber::Yarn yarn_tmp;

	const char* centerYarn1 = "centerYarn_ref.txt"; //TODO: maybe we can use z axis later
	yarn_tmp.yarnCenter(yarnfile1, centerYarn1);
	std::vector<yarnIntersect2D> pnts_ref;
	CrossSection cs2(yarnfile1, centerYarn1, ply_num, n, 100, pnts_ref);

	yarn_tmp.yarnCenter(yarnfile2, curveFile);
	std::vector<yarnIntersect2D> pnts_trans;
	CrossSection cs(yarnfile2, curveFile, ply_num, n, 100, pnts_trans);

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


	//optimize the solution and regularizing
	//cs.optimizeEllipses(ellipses, all_theta_R, new_ellipses, new_theta_R);
	new_ellipses = ellipses;
	new_theta_R = all_theta_R;

	////non-periodic theta
	//std::vector<float> theta;
	//std::vector<float> theta_new;
	//for (int i = 0; i < new_ellipses.size(); ++i)
	//	theta.push_back(new_ellipses[i].angle);
	//nonPeriodicTheta(theta, theta_new);
	//nonPeriodicTheta(new_theta_R, new_theta_R);
	//for (int i = 0; i < new_ellipses.size(); ++i)
	//	new_ellipses[i].angle = theta_new[i];

	FILE *fout;
	if (fopen_s(&fout, compressFile_vrtx, "wt") == 0) {
		fprintf_s(fout, "%d \n", ellipses.size());
		for (int i = 0; i < ellipses.size(); ++i) {
			fprintf_s(fout, "%.6f %.6f %.6f %.6f \n", new_ellipses[i].longR, new_ellipses[i].shortR, new_ellipses[i].angle, new_theta_R[i]);
		}
		fclose(fout);
	}

	FILE *fout0;
	if (fopen_s(&fout0, compressFile_seg, "wt") == 0) {
		//fprintf_s(fout0, "%d \n", ellipses.size());
		for (int i = 0; i < ellipses.size() - 1; ++i) {
			fprintf_s(fout0, "%.6f %.6f %.6f %.6f \n", (new_ellipses[i].longR+ new_ellipses[i+1].longR)/2.f, 
				(new_ellipses[i].shortR+ new_ellipses[i+1].shortR) / 2.f, 
				(new_ellipses[i].angle+ new_ellipses[i+1].angle) / 2.f, 
				(new_theta_R[i]+ new_theta_R[i+1]) / 2.f );
		}
		fclose(fout0);
	}

	//constant fitting
	std::cout << "Compression parameters are successfully written to the file!\n";

	//for debug: visualization
	const char* L2File = "../data/L2.txt";
	const char* refFile = "../data/allCrossSection2D_ref.txt";
	const char* deformedRefFile = "../data/allCrossSection2D_deformedRef.txt";
	const char* deformedFile = "../data/allCrossSection2D_deformed.txt";
	const float trimPercent = 0.0;
	plotIntersections(pnts_ref, refFile, trimPercent);
	std::vector<yarnIntersect2D> ref_deformed;
	deformRef(pnts_ref, ref_deformed, new_ellipses, new_theta_R);
	plotIntersections(ref_deformed, deformedRefFile, trimPercent);
	plotIntersections(pnts_trans, deformedFile, trimPercent);

	std::vector<float> L2;
	L2norm(ref_deformed, pnts_trans, L2, L2File);
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
		theta_new[i] = nextTheta(theta_new[i-1], theta_new[i]);
	}
	float k = floor(0.25*(theta_new[0] + theta_new[m-1]) / pi + 0.5)*2.0*pi;
	for (int i = 0; i < m; ++i) {
		theta_new[i] -= k;
	}
}

void constFitting_compParam( const std::vector<Ellipse> &ellipses, const std::vector<float> &theta_R,
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