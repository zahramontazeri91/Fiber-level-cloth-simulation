#include <random>
#include "../hermiteSeg.h"
#include "../hermiteCurve.h"
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
        curve.init(pts, norms, 10);

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
