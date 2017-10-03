#include "../hermiteSeg.h"
#include "hermiteTests.h"


void hermiteTest1()
{
    Eigen::Vector3d p0(0.0, 0.0, 0.0), p1(1.0, 1.0, 0.0);
    Eigen::Vector3d m0(0.0, 2.0, 0.0), m1(2.0, 0.0, 0.0);
    HermiteSpline hermite(p0, p1, m0, m1);

#if 1
    hermite.build(5);

    const int n = 100;
    Eigen::Vector3d vtx[n], tang[n], norm[n];
    hermite.output(n, vtx, tang, norm);
    FILE *fout;
    if ( fopen_s(&fout, "junk.txt", "wt") == 0 ) {
        for ( int i = 0; i < n; ++i ) {
            fprintf_s(fout, "%.4lf %.4lf %.4lf ", vtx[i][0], vtx[i][1], vtx[i][2]);
            fprintf_s(fout, "%.4lf %.4lf %.4lf ", tang[i][0], tang[i][1], tang[i][2]);
            fprintf_s(fout, "%.4lf %.4lf %.4lf\n", norm[i][0], norm[i][1], norm[i][2]);
        }
        fclose(fout);
    }
#else
    hermite.build(100);

    double t = hermite.arcLengthInvApprox(0.567);
    printf("%.6lf\n", hermite.arcLengthApprox(t));
#endif
}
