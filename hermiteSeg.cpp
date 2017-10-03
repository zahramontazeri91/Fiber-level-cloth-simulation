//#include "stdafx.h"
#include <algorithm>
#include "hermiteSeg.h"

#define HERMITE_EPS     1e-5
#define POLYSOLVER_EPS  1e-8


template <int n>
struct polySolver
{
    typedef typename Eigen::Matrix<double, 1, n + 1> InputType;

    static int run(const InputType &c, double *roots)
    {
        const Eigen::Matrix<double, 1, n> &c1 = c.block<1, n>(0, 1);
        if ( std::abs(c[0]) < POLYSOLVER_EPS ) {
            //SLog(EWarn, "Unstable polynomial solve: c5 = %.2e\n", c[0]);
            return polySolver<n - 1>::run(c1, roots);
        }

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


double HermiteSpline::arcLengthApprox(double t, int subdiv) const
{
    double res = 0.0;
    Eigen::Vector3d last = p0, cur;
    double t0 = 0.0;
    for ( int i = 1; i <= subdiv; ++i ) {
        double t1 = static_cast<double>(i)/subdiv;
        cur = eval(t1);
        if ( t > t1 ) {
            res += (last - cur).norm();
            last = cur; t0 = t1;
        }
        else {
            res += (last - cur).norm()*(t - t0)/(t1 - t0);
            break;
        }
    }
    return res;
}


double HermiteSpline::arcLengthInvApprox(double len, int subdiv) const
{
	double L = 0.0, R = 1.0;
	while (R - L > HERMITE_EPS) {
		double mid = 0.5*(L + R);
		double val = arcLengthApprox(mid, subdiv);
		if (val > len)
			R = mid;
		else
			L = mid;
	}
	return 0.5*(L + R);
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


void HermiteSpline::output(int n, Eigen::Vector3d *buffer) const
{
	assert(n > 1);
	for (int i = 0; i < n; ++i) {
		double t = static_cast<double>(i) / (n - 1);
		buffer[i] = eval(t);
	}
}

void HermiteSpline_multiSeg::init_splines(const char* filename) {
	// read the spline parameters file
	std::ifstream fin;
	if (filename != NULL)
		fin.open(filename);

	std::string line;
	std::getline(fin, line);
	m_spline_seg = atof(line.c_str());

	for (int i = 0; i < m_spline_seg; i++) {

		std::getline(fin, line);
		std::vector<std::string> splits = split(line, ' ');
		const Eigen::Vector3d p0(atof(splits[0].c_str()), atof(splits[1].c_str()), atof(splits[2].c_str()));

		std::getline(fin, line);
		splits = split(line, ' ');
		const Eigen::Vector3d p1(atof(splits[0].c_str()), atof(splits[1].c_str()), atof(splits[2].c_str()));

		std::getline(fin, line);
		splits = split(line, ' ');
		const Eigen::Vector3d m0(atof(splits[0].c_str()), atof(splits[1].c_str()), atof(splits[2].c_str()));

		std::getline(fin, line);
		splits = split(line, ' ');
		const Eigen::Vector3d m1(atof(splits[0].c_str()), atof(splits[1].c_str()), atof(splits[2].c_str()));

		// initialize the splines
		HermiteSpline spline(p0, p1, m0, m1);
		m_splines.push_back(spline);
	}
}

Eigen::Vector3d HermiteSpline_multiSeg::eval(double t) {
	int spline_id = int(t); // curve from t = 0 to t = 1 defined by first segment and so on
	if (t == m_spline_seg)
		return m_splines[spline_id - 1].eval(1.0);
	return m_splines[spline_id].eval(t - double(spline_id));
}

Eigen::Vector3d HermiteSpline_multiSeg::evalTangent(double t) {
	int spline_id = int(t); // curve from t = 0 to t = 1 defined by first segment and so on
	if (t == m_spline_seg)
		return m_splines[spline_id - 1].evalTangent(1.0);
	return m_splines[spline_id].evalTangent(t - double(spline_id));
}

Eigen::Vector3d HermiteSpline_multiSeg::evalCurvature(double t) {
	int spline_id = int(t);
	if (t == m_spline_seg)
		return m_splines[spline_id - 1].evalCurvature(1.0);
	return m_splines[spline_id].evalCurvature(t - double(spline_id));
}

double HermiteSpline_multiSeg::arcLengthInvApprox(double len, int subdiv) {

	int spline_id = int(len);
	if (len == m_spline_seg)
		return m_splines[spline_id - 1].arcLengthInvApprox(1.0, subdiv);
	return m_splines[spline_id].arcLengthInvApprox(len - double(spline_id), subdiv);
}

std::vector<double> HermiteSpline_multiSeg::segLengths(int subdiv) {
	std::vector<double> length_list;
	for (int i = 0; i < m_spline_seg; i++) {
		double length_seg_i = m_splines[i].totalLength(subdiv);
		length_list.push_back(length_seg_i);
	}
	return length_list;
}

double HermiteSpline_multiSeg::totalLength(int subdiv) {
	std::vector<double> length_list = segLengths(subdiv);
	double sum = 0.0;
	for (int i = 0; i < length_list.size(); i++) {
		sum += length_list[i];
	}
	return sum;
}

int HermiteSpline_multiSeg::findSegId(double curve_length, int subdiv) { //TODO:move length_list to class member
	int segId = 0;
	double length_sofar = 0.0;
	std::vector<double> length_list = segLengths(subdiv);
	
	for (int i = 0; i < length_list.size(); i++ ){
		if ((curve_length - length_sofar) >= 0) {
			length_sofar += length_list[i];
			segId = i ;
			continue;
		}
		else {
			segId = i - 1;
			break;
		}
	}
	return segId;
}