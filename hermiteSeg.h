#pragma once
#ifndef __HERMITE_SPLINE_H_
#define __HERMITE_SPLINE_H_

#include <vector>
#include <sstream>
#include <Eigen/Dense>

#define HERMITE_STRICT_PROJECTION
//#define HERMITE_ENABLE_BRUTEFORCE

#define HERMITE_EPS     1e-5


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

    void build(int subdiv, Eigen::Vector3d norm0 = Eigen::Vector3d::Zero(),
        Eigen::Vector3d norm1 = Eigen::Vector3d::Zero());

	/* get the position for some given 0 <= t <= 1*/
    inline Eigen::Vector3d eval(double t) const
    {
        return ((u1*t + u2)*t + m0)*t + p0; // the hermite formula
    }

	/* get the tangent for some given 0 <= t <= 1 */
    inline Eigen::Vector3d evalTangent(double t, bool normalize = true) const
    {
        Eigen::Vector3d ret = (u1*3.0*t + u2*2.0)*t + m0;
        if ( normalize ) {
            //assert(ret.norm() > HERMITE_EPS); 
            ret.normalize();
        }
        return ret;
    }

	/* get the curvature (i.e., derivative of the velocity) at t */
	inline Eigen::Vector3d evalCurvature(double t) const
	{
		Eigen::Vector3d ret = u1*6.0*t + u2*2.0;
		return ret;
	}
	
	inline Eigen::Vector3d rotateTang(const Eigen::Vector3d &v) const {
		Eigen::Matrix3d Rx, Ry, Rz;
		Rx << 1, 0, 0,
			0, 0, -1,
			0, 1, 0;
		Ry << 0, 0, 1,
			0, 1, 0,
			-1, 0, 0;
		Rz << 0, -1, 0,
			1, 0, 0,
			0, 0, 1;
		return (Rz*v == v ? Ry*v : Rz*v);
	}

    /* get the principle normal at t (new) */ 
    inline Eigen::Vector3d evalPrincipalNormal(double t, bool normalize = true) const
    {

		Eigen::Vector3d ret;
        Eigen::Vector3d q = evalCurvature(t), v = evalTangent(t);

		q.normalize();

		if (normalize) {
			assert(q.norm() > HERMITE_EPS);
			assert(v.norm() > HERMITE_EPS);

			q.normalize(); //normalize the curvature before VxQxV 
			v.normalize();
		}

        ret = v.cross(q).cross(v);
        if ( normalize ) {
			assert(ret.norm() > HERMITE_EPS && "Either normal or tangent is zero!");
			ret.normalize();            
        }
		assert(std::abs(ret.dot(v)) < HERMITE_EPS && "normal and tangent are not perpendicular!");
        return ret;
    }


	/* Get the total length of the spline */
	inline double totalLength() {
		return arcLengthApprox(1.0);
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

    Eigen::Vector3d evalNormal(double t) const;

    bool project(const Eigen::Vector3d &x, double &t, double &dist, Eigen::Vector3d *q = NULL) const;
#ifdef HERMITE_ENABLE_BRUTEFORCE
    double project_BruteForce(const Eigen::Vector3d &x, double &dist, Eigen::Vector3d *q = NULL) const;
#endif

    double errorBound() const;
#ifdef HERMITE_ENABLE_BRUTEFORCE
    double errorBound_BruteForce() const;
#endif

	/* get the arc length at t */
    double arcLengthApprox(double t) const;

	/* get the t value that gives the arc length len */
    double arcLengthInvApprox(double len) const;

    double subdivideN(int n, HermiteSpline *splines, double *error = NULL) const;
    double subdivideA(double maxError, std::vector<HermiteSpline> &results) const;

	// return the normal for the other end given tangents for both ends
    static Eigen::Vector3d computeRotatedNormal(const Eigen::Vector3d &tang0, const Eigen::Vector3d &tang1,
        const Eigen::Vector3d norm0);
	static Eigen::Vector3d computeRotatedNormal_backward(const Eigen::Vector3d &tang0, const Eigen::Vector3d &tang1,
		const Eigen::Vector3d norm1);

    void output(int n, Eigen::Vector3d *bufferPosition,
        Eigen::Vector3d *bufferTangent = NULL, Eigen::Vector3d *bufferNormal = NULL) const;

	void getNorms(std::vector<Eigen::Vector3d> &m_norms) {
		m_norms = norms;
	}
protected:
    double subdivideAInternal(double maxError, std::vector<HermiteSpline> &results) const;

    Eigen::Vector3d p0, p1, m0, m1;
    Eigen::Vector3d u1, u2;

    int subdiv;
    std::vector<double> lens;
    std::vector<Eigen::Vector3d> tangents, norms;
};

#endif
