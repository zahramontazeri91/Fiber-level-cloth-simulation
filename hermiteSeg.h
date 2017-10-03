#pragma once
#ifndef __HERMITE_SPLINE_H_
#define __HERMITE_SPLINE_H_

#include <vector>
#include <sstream>
#include <Eigen/Dense>
#include "Util.h"

#define HERMITE_STRICT_PROJECTION
//#define HERMITE_ENABLE_BRUTEFORCE


class HermiteSpline
{
public:
	/* Single segment of the spline given by the position (p0, p1) and tangent (m0, m1) of the two endpoints */
    HermiteSpline(const Eigen::Vector3d &_p0, const Eigen::Vector3d &_p1,
        const Eigen::Vector3d &_m0, const Eigen::Vector3d &_m1)
    {
        init(_p0, _p1, _m0, _m1);
    }

	/* creating a Hermite spline segment */
    void init(const Eigen::Vector3d &_p0, const Eigen::Vector3d &_p1,
        const Eigen::Vector3d &_m0, const Eigen::Vector3d &_m1)
    {
        p0 = _p0; p1 = _p1; m0 = _m0; m1 = _m1;
        u1 = (p0 - p1)*2.0 + m0 + m1;
        u2 = (p0 - p1)*(-3.0) - m0*2.0 - m1;
    }

    inline HermiteSpline()
    {
        init(Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(1.0, 0.0, 0.0),
             Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(0.0, -1.0, 0.0));
    }

    inline HermiteSpline(const HermiteSpline &spline)
    {
        init(spline.p0, spline.p1, spline.m0, spline.m1);
    }

	/* get the position for some given 0 <= t <= 1*/
    inline Eigen::Vector3d eval(double t) const
    {
        return ((u1*t + u2)*t + m0)*t + p0; //the hermite formula
    }

	/* get the tangent for some given 0 <= t <= 1 */
    inline Eigen::Vector3d evalTangent(double t) const
    {
        return (u1*3.0*t + u2*2.0)*t + m0;
    }

	/* get the derivative of the velocity */
	inline Eigen::Vector3d evalCurvature(double t) const
	{
		return (u1*6.0*t + u2*2.0);
	}

	/* Get the total length of the spline */
	inline float totalLength(int subdiv) {
		return arcLengthApprox(1, subdiv);
	}

    inline std::string toString() const
    {
        std::ostringstream oss;
        oss.precision(4);
        oss.setf(std::ios::scientific);

        oss << p0[0] << ' ' << p0[1] << ' ' << p0[2] << ' '
            << p1[0] << ' ' << p1[1] << ' ' << p1[2] << ' '
            << m0[0] << ' ' << m0[1] << ' ' << m0[2] << ' '
            << m1[0] << ' ' << m1[1] << ' ' << m1[2];
        return oss.str();
    }


    bool project(const Eigen::Vector3d &x, double &t, double &dist, Eigen::Vector3d *q = NULL) const;
#ifdef HERMITE_ENABLE_BRUTEFORCE
    double project_BruteForce(const Eigen::Vector3d &x, double &dist, Eigen::Vector3d *q = NULL) const;
#endif

    double errorBound() const;
#ifdef HERMITE_ENABLE_BRUTEFORCE
    double errorBound_BruteForce() const;
#endif

	/* get the arc length at t */
    double arcLengthApprox(double t, int subdiv) const;

	/* get the t value that gives the arc length len*/

    double arcLengthInvApprox(double len, int subdiv) const; // Slow implementation

    double subdivideN(int n, HermiteSpline *splines, double *error = NULL) const;
    double subdivideA(double maxError, std::vector<HermiteSpline> &results) const;

    void output(int n, Eigen::Vector3d *buffer) const;

	/* return z value for the starting point */
	inline double get_start_z() {
		return p0[2];
	}
protected:
    double subdivideAInternal(double maxError, std::vector<HermiteSpline> &results) const;

    Eigen::Vector3d p0, p1, m0, m1;
    Eigen::Vector3d u1, u2;
};

#endif
