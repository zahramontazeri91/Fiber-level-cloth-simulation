#include <random>
#include "../hermiteSeg.h"
#include "../hermiteCurve.h"
#include "../Fiber.h"
#include "../fitting.h"
#include "../crossSection.h"
#include "hermiteTests.h"


void hermiteTest1()
{
    // Testing single segment
    {
        Eigen::Vector3d p0(0.0, 0.0, 0.0), p1(1.0, 1.0, 0.0);
        Eigen::Vector3d m0(0.0, 2.0, 0.0), m1(2.0, 0.0, 0.0);
        HermiteSpline hermite(p0, p1, m0, m1);
        hermite.build(5);

        const int n = 100;
        Eigen::Vector3d vtx[n], tang[n], norm[n];
        hermite.output(n, vtx, tang, norm);
        FILE *fout;
        if ( fopen_s(&fout, "junk_single.txt", "wt") == 0 ) {
            for ( int i = 0; i < n; ++i ) {
                fprintf_s(fout, "%.4lf %.4lf %.4lf ", vtx[i][0], vtx[i][1], vtx[i][2]);
                fprintf_s(fout, "%.4lf %.4lf %.4lf ", tang[i][0], tang[i][1], tang[i][2]);
                fprintf_s(fout, "%.4lf %.4lf %.4lf\n", norm[i][0], norm[i][1], norm[i][2]);
            }
            fclose(fout);
        }
    }

    // Testing entire curve
    {
        const int nSeg = 10;
        std::vector<Eigen::Vector3d> pts;
        for ( int i = 0; i <= nSeg; ++i )
            pts.push_back(Eigen::Vector3d(static_cast<double>(i), static_cast<double>(i % 2), std::sin(static_cast<double>(i))));
        HermiteCurve curve;
        curve.init(pts, 10);

        const int n = 200;
        Eigen::Vector3d vtx[n], tang[n], norm[n];
        curve.output(n, vtx, tang, norm);
        FILE *fout;
        if ( fopen_s(&fout, "junk_multiple.txt", "wt") == 0 ) {
            for ( int i = 0; i < n; ++i ) {
                fprintf_s(fout, "%.4lf %.4lf %.4lf ", vtx[i][0], vtx[i][1], vtx[i][2]);
                fprintf_s(fout, "%.4lf %.4lf %.4lf ", tang[i][0], tang[i][1], tang[i][2]);
                fprintf_s(fout, "%.4lf %.4lf %.4lf\n", norm[i][0], norm[i][1], norm[i][2]);
            }
            fclose(fout);
        }
    }
}


void hermiteTest2()
{
    // Testing single segment
    {
        Eigen::Vector3d p0(0.0, 0.0, 0.0), p1(1.0, 1.0, 0.0);
        Eigen::Vector3d m0(0.0, 2.0, 0.0), m1(2.0, 0.0, 0.0);
        HermiteSpline hermite(p0, p1, m0, m1);

        hermite.build(100);

        double t = hermite.arcLengthInvApprox(0.567);
        printf("%.6lf\n", hermite.arcLengthApprox(t));
    }

    // Testing entire curve
    {
        const int nSeg = 40;

        std::vector<Eigen::Vector3d> pts;
        for ( int i = 0; i <= nSeg; ++i )
            pts.push_back(Eigen::Vector3d(static_cast<double>(i), static_cast<double>(i % 2), std::sin(static_cast<double>(i))));
        HermiteCurve curve;
        curve.init(pts, 10);

        std::uniform_real_distribution<> distrb;
        std::mt19937_64 engine;

        for ( int i = 0; i < 100; ++i ) {
            double t0 = distrb(engine)*static_cast<double>(nSeg);
            double len = curve.arcLengthApprox(t0);
            double t1 = curve.arcLengthInvApprox(len);
            assert(std::abs(t0 - t1) < HERMITE_EPS);
        }
    }
}


void hermiteTest3()
{
    // Testing entire curve
    {
        const int nSeg = 50;
        std::vector<Eigen::Vector3d> pts, norms;
        for ( int i = 0; i <= nSeg; ++i ) {
            pts.push_back(Eigen::Vector3d(0.0, 0.0, 0.1*i));

            double angle = 4.0*std::acos(0.0)*i/10.0;
            norms.push_back(Eigen::Vector3d(std::cos(angle), std::sin(angle), 0.0));
        }
        HermiteCurve curve;
        curve.init_norm(pts, norms, 10);

        const int n = 200;
        Eigen::Vector3d vtx[n], tang[n], norm[n];
        curve.output(n, vtx, tang, norm);
        FILE *fout;
        if ( fopen_s(&fout, "junk_multiple_2.txt", "wt") == 0 ) {
            for ( int i = 0; i < n; ++i ) {
                fprintf_s(fout, "%.4lf %.4lf %.4lf ", vtx[i][0], vtx[i][1], vtx[i][2]);
                fprintf_s(fout, "%.4lf %.4lf %.4lf ", tang[i][0], tang[i][1], tang[i][2]);
                fprintf_s(fout, "%.4lf %.4lf %.4lf\n", norm[i][0], norm[i][1], norm[i][2]);
            }
            fclose(fout);
        }
    }
}


void vis_allFrames() {
//	Fiber::Yarn yarn;
//	const char* configfile = "config.txt";
//	yarn.parse(configfile);
//	const int n = yarn.getStepNum();
//	const int ply_num = yarn.getPlyNum();
//
//	for (int i = 10; i < 30; ++i) {
//		std::string s;
//		if (i < 10)
//			s = "D:/sandbox/fiberSimulation/hairs/scaledHairs/frame0000" + std::to_string(i) + "_scaled.txt";
//		else
//			s = "D:/sandbox/fiberSimulation/hairs/scaledHairs/frame000" + std::to_string(i) + "_scaled.txt";
//		const char* yarnfile1 = s.c_str();
//		std::ifstream fin(yarnfile1);
//		std::cout << yarnfile1 << std::endl;
//		assert(fin.is_open() && "file wasn't found!\n");
//
//		Fiber::Yarn yarn_tmp;
//		const char* centerYarn1 = "centerYarn_ref.txt";
//		yarn_tmp.yarnCenter(yarnfile1, centerYarn1);
//		std::vector<yarnIntersect2D> pnts;
//		CrossSection cs2(yarnfile1, centerYarn1, ply_num, n, 100, pnts);
//
//		std::string t = "../data/allCrossSection2D_deformed_frame" + std::to_string(i) + ".txt";
//		const char* deformed = t.c_str();
//		plotIntersections(pnts, deformed, 0.3);
//	}
}