
#include <algorithm>
#include "hermiteSeg.h"

#include <iostream>

#define POLYSOLVER_EPS  1e-8


template <int n>
struct polySolver
{
    typedef typename Eigen::Matrix<double, 1, n + 1> InputType;

    static int run(const InputType &c, double *roots)
    {
        const Eigen::Matrix<double, 1, n> &c1 = c.block<1, n>(0, 1);
        if ( std::abs(c[0]) < POLYSOLVER_EPS ) return polySolver<n - 1>::run(c1, roots);

        Eigen::Matrix<double, n, n> A = Eigen::Matrix<double, n, n>::Zero();
        A.diagonal(-1) = Eigen::Matrix<double, n - 1, 1>::Ones();
        A.row(0) = -c1/c[0];

        const auto &eigVals = A.eigenvalues();
        int nroots = 0;
        for ( int i = 0; i < eigVals.rows(); ++i )
            if ( std::abs(eigVals(i, 0).imag()) < POLYSOLVER_EPS )
                roots[nroots++] = eigVals(i, 0).real();
        if ( nroots ) std::sort(roots, roots + nroots);
        return nroots;
    }
};


template <>
struct polySolver<1>
{
    typedef Eigen::Matrix<double, 1, 2> InputType;

    static int run(const InputType &c, double *roots)
    {
        int nroots = 0;
        if ( std::abs(c[0]) > POLYSOLVER_EPS ) {
            roots[0] = -c[1]/c[0];  //solve ax+b = 0
            nroots = 1;
        }
        return nroots;
    }
};

#if 0
void HermiteSpline::build(int _subdiv, Eigen::Vector3d _norm0, Eigen::Vector3d _norm1)
{
	subdiv = _subdiv;

	std::vector<Eigen::Vector3d> positions(subdiv + 1);
	tangents.resize(subdiv + 1);
	for (int i = 0; i <= subdiv; ++i) {
		double t = static_cast<double>(i) / subdiv;
		positions[i] = eval(t);
		tangents[i] = evalTangent(t);
	}

	lens.resize(subdiv + 1);
	lens[0] = 0.0;
	for (int i = 1; i <= subdiv; ++i)
		lens[i] = lens[i - 1] + (positions[i - 1] - positions[i]).norm();

	norms.resize(subdiv + 1);
	norms[0] = _norm0.norm() > HERMITE_EPS ? _norm0.normalized() : evalPrincipalNormal(0.0);
	if (_norm1.norm() > HERMITE_EPS) {
		norms[subdiv] = _norm1.normalized();
		for (int i = 1; i < subdiv; ++i) {
			double t = static_cast<double>(i) / subdiv;
			norms[i] = ((1.0 - t)*_norm0 + t*_norm1).normalized();
		}
	}
	else {
		for (int i = 1; i <= subdiv; ++i)
			norms[i] = computeRotatedNormal(tangents[i - 1], tangents[i], norms[i - 1]);
	}
}
#endif

#if 1
void HermiteSpline::build(int _subdiv, Eigen::Vector3d _norm0, Eigen::Vector3d _norm1)
{
    subdiv = _subdiv;

    std::vector<Eigen::Vector3d> positions(subdiv + 1);
    tangents.resize(subdiv + 1);
    for ( int i = 0; i <= subdiv; ++i ) {
        double t = static_cast<double>(i)/subdiv;
        positions[i] = eval(t);
        tangents[i] = evalTangent(t);
    }

    lens.resize(subdiv + 1);
    lens[0] = 0.0;
    for ( int i = 1; i <= subdiv; ++i )
        lens[i] = lens[i - 1] + (positions[i - 1] - positions[i]).norm();

    norms.resize(subdiv + 1);

	//interpolate for specific segment for the case of excisting normals for both ends
    if (_norm0.norm() >= HERMITE_EPS && _norm1.norm() >= HERMITE_EPS ) {
		norms[0] = _norm0.normalized();
        norms[subdiv] = _norm1.normalized();
        for ( int i = 1; i < subdiv; ++i ) {
            double t = static_cast<double>(i)/subdiv;
            norms[i] = ((1.0 - t)*_norm0 + t*_norm1).normalized();
        }
    }
	//use forward-RM if the first normal is given 
    else if (_norm0.norm() > HERMITE_EPS && _norm1.norm() < HERMITE_EPS) {
		norms[0] = _norm0.normalized();
		for (int i = 1; i <= subdiv; ++i) {
			norms[i] = computeRotatedNormal(tangents[i - 1], tangents[i], norms[i - 1]);
		}
    }
	//use backward-RM if the first normal is given 
	else {
		norms[subdiv] = _norm1.normalized();
		for (int i = subdiv - 1; i <= 0; --i) {
			norms[i] = computeRotatedNormal_backward(tangents[i], tangents[i + 1], norms[i + 1]);
			assert(std::abs(norms[i].norm() - 1.0) < HERMITE_EPS);
		}
	}

	for (int i = 0; i <= subdiv; ++i) {
		norms[i].normalize();
	}
}
#endif

Eigen::Vector3d HermiteSpline::evalNormal(double t) const
{
	
    if ( subdiv ) {

		if ( t < HERMITE_EPS ) return norms[0].normalized();
        if ( t > 1.0 - HERMITE_EPS ) return norms[subdiv].normalized();
        int idx = static_cast<int>(std::floor(t*subdiv));
        Eigen::Vector3d tang0 = evalTangent(t);
        Eigen::Vector3d norm0 = computeRotatedNormal(tangents[idx], tang0, norms[idx]),
                        norm1 = computeRotatedNormal(tangents[idx + 1], tang0, norms[idx + 1]);
        double w = (t*subdiv) - idx;
        Eigen::Vector3d ret = (1.0 - w)*norm0 + norm1;
        ret = tang0.cross(ret).cross(tang0);
        assert(ret.norm() > HERMITE_EPS);
        ret.normalize();
        assert(std::abs(ret.dot(tang0)) < HERMITE_EPS);
		assert(std::abs(ret.norm() - 1.0) < HERMITE_EPS);
        return ret;
    }
    else {
        fprintf(stderr, "Error: normal uninitialized!\n");
        return Eigen::Vector3d::Zero();
    }
}


bool HermiteSpline::project(const Eigen::Vector3d &x, double &t, double &dist, Eigen::Vector3d *q) const
{
    const Eigen::Vector3d e = p0 - x;
    Eigen::Matrix<double, 1, 6> c;
    c << 3.0*u1.squaredNorm(), 5.0*u1.dot(u2), 4.0*u1.dot(m0) + 2.0*u2.squaredNorm(),
         3.0*u1.dot(e) + 3.0*u2.dot(m0), 2.0*u2.dot(e) + m0.squaredNorm(), e.dot(m0);
    double roots[5];
    int nroots = polySolver<5>::run(c, roots);

    double d0, d1;
#ifdef HERMITE_STRICT_PROJECTION
    bool done = false;
    dist = std::numeric_limits<double>::infinity();
#else
    d0 = distance(p0, x); d1 = distance(p1, x);
    if ( d0 < d1 ) {
        dist = d0; t = 0.0;
        if ( q ) *q = p0;
    }
    else {
        dist = d1; t = 1.0;
        if ( q ) *q = p1;
    }
#endif

    for ( int i = 0; i < nroots; ++i )
#ifdef HERMITE_STRICT_PROJECTION
        if ( roots[i] > -HERMITE_EPS && roots[i] < 1.0 + HERMITE_EPS ) {
#else
        if ( roots[i] > HERMITE_EPS && roots[i] < 1.0 - HERMITE_EPS ) {
#endif
            d1 = roots[i];
            Eigen::Vector3d q0 = eval(d1);
            d0 = (q0 - x).norm();
            if ( d0 < dist ) {
                dist = d0; t = d1;
                if ( q ) *q = q0;
#ifdef HERMITE_STRICT_PROJECTION
                done = true;
#endif
            }
        }

#ifdef HERMITE_STRICT_PROJECTION
    return done;
#else
    return true;
#endif
}


#ifdef HERMITE_ENABLE_BRUTEFORCE
double HermiteSpline::project_BruteForce(const Point3d &x, double &dist, Point3d *q) const
{
    const int N = 100000;
    dist = distance(p0, x);
    double ans = 0.0;
    for ( int i = 1; i <= N; ++i ) {
        double t = static_cast<double>(i)/N;
        Point3d q0 = eval(t);
        double d0 = distance(q0, x);
        if ( d0 < dist ) {
            dist = d0; ans = t;
            if ( q ) *q = q0;
        }
    }
    return ans;
}
#endif


double HermiteSpline::errorBound() const
{
    Eigen::Vector3d u3 = m0 - (p1 - p0);
    Eigen::Matrix<double, 1, 6> c;
    c << 3.0*u1.squaredNorm(), 5.0*u1.dot(u2), 4.0*u1.dot(u3) + 2.0*u2.squaredNorm(),
        3.0*u2.dot(u3), u3.squaredNorm(), 0.0;
    double roots[5];
    int nroots = polySolver<5>::run(c, roots);

    double best = 0.0;
    for ( int i = 0; i < nroots; ++i )
        if ( roots[i] > HERMITE_EPS && roots[i] < 1.0 - HERMITE_EPS ) {
            double t = roots[i];
            double dist = (eval(t) - (p0 + (p1 - p0)*t)).norm();
            if ( dist > best ) best = dist;
        }
    return best;
}


#ifdef HERMITE_ENABLE_BRUTEFORCE
double HermiteSpline::errorBound_BruteForce() const
{
    const int N = 100000;
    double best = 0.0;
    for ( int i = 1; i < N; ++i ) {
        double t = static_cast<double>(i)/N;
        double d0 = distance(eval(t), p0 + (p1 - p0)*t);
        if ( d0 > best ) best = d0;
    }
    return best;
}
#endif


double HermiteSpline::arcLengthApprox(double t) const
{
    if ( t < HERMITE_EPS ) return 0.0;
    if ( t > 1.0 - HERMITE_EPS ) return lens[subdiv];

    int i = static_cast<int>(std::floor(t*subdiv));
    double w = t*subdiv - i;
    return lens[i] + w*(lens[i + 1] - lens[i]);
}


double HermiteSpline::arcLengthInvApprox(double len) const
{
    if ( len < HERMITE_EPS ) return 0.0;
    if ( len > lens[subdiv] - HERMITE_EPS ) return 1.0;

    int i = static_cast<int>(std::lower_bound(lens.begin(), lens.end(), len) - lens.begin());
    return (i - 1 + (len - lens[i - 1])/(lens[i] - lens[i - 1]))/subdiv;
}


double HermiteSpline::subdivideN(int n, HermiteSpline *splines, double *error) const
{
	double maxError = 0.0;
	for (int i = 0; i < n; ++i) {
		double t0 = static_cast<double>(i) / n,
			t1 = static_cast<double>(i + 1) / n;
		splines[i].init(eval(t0), eval(t1),
			evalTangent(t0) / static_cast<double>(n), evalTangent(t1) / static_cast<double>(n));
		double errVal = splines[i].errorBound();
		maxError = std::max(maxError, errVal);
		if (error) error[i] = errVal;
	}
	return maxError;
}


double HermiteSpline::subdivideA(double maxError, std::vector<HermiteSpline> &results) const
{
	results.clear();
	double errVal = errorBound();
	if (errVal < maxError)
		results.push_back(*this);
	else
		errVal = subdivideAInternal(maxError, results);
	return errVal;
}


double HermiteSpline::subdivideAInternal(double maxError, std::vector<HermiteSpline> &results) const
{
	HermiteSpline splines[2];
	double errVal[2], ret = 0.0;

	subdivideN(2, splines, errVal);
	for (int i = 0; i < 2; ++i)
		if (errVal[i] < maxError) {
			results.push_back(splines[i]);
			ret = std::max(ret, errVal[i]);
		}
		else
			ret = std::max(ret, splines[i].subdivideAInternal(maxError, results));
			return ret;
}


void HermiteSpline::output(int n, Eigen::Vector3d *bufferPosition, Eigen::Vector3d *bufferTangent, Eigen::Vector3d *bufferNormal) const
{
    assert(n > 1);
    for ( int i = 0; i < n; ++i ) {
        double t = static_cast<double>(i)/(n - 1);
        bufferPosition[i] = eval(t);
        if ( bufferTangent ) bufferTangent[i] = evalTangent(t);
        if ( bufferNormal ) bufferNormal[i] = evalNormal(t);
    }
}


Eigen::Vector3d HermiteSpline::computeRotatedNormal(const Eigen::Vector3d &tang0, const Eigen::Vector3d &tang1, const Eigen::Vector3d norm0)
{
	//std::cout << std::abs(norm0.norm() - 1.0) << " ******** \n " << norm0 << std::endl;

    assert(std::abs(tang0.norm() - 1.0) < HERMITE_EPS);
    assert(std::abs(tang1.norm() - 1.0) < HERMITE_EPS);
	assert(std::abs(norm0.norm() - 1.0) < HERMITE_EPS);

    double val = tang0.dot(tang1);

	//make sure it doesn't always return without rotating the normal
	if (val > 1.0 - 1e-7)
		return norm0;

    else if ( val < -1.0 + 1e-7) {
        fprintf(stderr, "Warning: oppositely pointing tangents\n");
        return -norm0;
    }

    Eigen::Matrix3d m = Eigen::AngleAxisd(std::acos(val), tang0.cross(tang1).normalized()).toRotationMatrix();
    assert((m*tang0 - tang1).norm() < HERMITE_EPS);
    return (m*norm0).normalized();
}

Eigen::Vector3d HermiteSpline::computeRotatedNormal_backward(const Eigen::Vector3d &tang0, const Eigen::Vector3d &tang1, const Eigen::Vector3d norm1)
{
	assert(std::abs(tang0.norm() - 1.0) < HERMITE_EPS);
	assert(std::abs(tang1.norm() - 1.0) < HERMITE_EPS);
	assert(std::abs(norm1.norm() - 1.0) < HERMITE_EPS);

	double val = tang0.dot(tang1);
	if (val > 1.0 - HERMITE_EPS)
		return norm1;
	else if (val < -1.0 + HERMITE_EPS) {
		fprintf(stderr, "Warning: oppositely pointing tangents\n");
		return -norm1;
	}

	Eigen::Matrix3d m = Eigen::AngleAxisd(std::acos(val), tang0.cross(tang1).normalized()).toRotationMatrix();
	assert((m*tang0 - tang1).norm() < HERMITE_EPS);
	assert(std::abs((m*norm1).norm() - 1.0) < HERMITE_EPS);
	return (m*norm1).normalized();
}