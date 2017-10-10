#include <algorithm>
#include "hermiteCurve.h"


void HermiteCurve::init(const char* filename, int subdiv) {
    assert(filename);
	std::ifstream fin(filename);
    assert(fin.is_open());

    int n;
    fin >> n;

    std::vector<Eigen::Vector3d> pts(n);
    for ( int i = 0; i < n; ++i ) fin >> pts[i][0] >> pts[i][1] >> pts[i][2];

    init(pts, subdiv);
}


void HermiteCurve::init(const std::vector<Eigen::Vector3d> &pts, int subdiv)
{
    assert(pts.size()>2);

    m_spline_seg = static_cast<int>(pts.size()) - 1;
    m_splines.resize(m_spline_seg);
    for ( int i = 0; i < m_spline_seg; ++i ) {
        Eigen::Vector3d m0, m1;
        if ( i == 0 )
            m0 = pts[i + 1] - pts[i];
        else
            m0 = (pts[i + 1] - pts[i - 1])*0.5;
        if ( i == m_spline_seg - 1 )
            m1 = pts[i + 1] - pts[i];
        else
            m1 = (pts[i + 2] - pts[i])*0.5;
        m_splines[i].init(pts[i], pts[i + 1], m0, m1);
    }

    m_splines[0].build(subdiv, m_splines[0].evalPrincipalNormal(0.0));

    for ( int i = 1; i < m_spline_seg; ++i )
        m_splines[i].build(subdiv, m_splines[i - 1].evalNormal(1.0));
	
    m_lens.resize(m_spline_seg);
    for ( int i = 0; i < m_spline_seg; ++i ) {
        m_lens[i] = m_splines[i].totalLength();
        if ( i ) m_lens[i] += m_lens[i - 1];
    }
}


Eigen::Vector3d HermiteCurve::eval(double t) const
{
    assert(m_spline_seg > 0);
    if ( t < HERMITE_EPS ) return m_splines[0].eval(0.0);
    if ( t > m_spline_seg - HERMITE_EPS ) return m_splines[m_spline_seg - 1].eval(1.0);
    return m_splines[static_cast<int>(t)].eval(t - std::floor(t));
}


Eigen::Vector3d HermiteCurve::evalTangent(double t, bool normalize) const
{
    assert(m_spline_seg > 0);
    if ( t < HERMITE_EPS ) return m_splines[0].evalTangent(0.0, normalize);
    if ( t > m_spline_seg - HERMITE_EPS ) return m_splines[m_spline_seg - 1].evalTangent(1.0, normalize);
    return m_splines[static_cast<int>(t)].evalTangent(t - std::floor(t), normalize);
}


Eigen::Vector3d HermiteCurve::evalCurvature(double t) const
{
    assert(m_spline_seg > 0);
    if ( t < HERMITE_EPS ) return m_splines[0].evalCurvature(0.0);
    if ( t > m_spline_seg - HERMITE_EPS ) return m_splines[m_spline_seg - 1].evalCurvature(1.0);
    return m_splines[static_cast<int>(t)].evalCurvature(t - std::floor(t));
}


Eigen::Vector3d HermiteCurve::evalNormal(double t) const
{
    assert(m_spline_seg > 0);
    if ( t < HERMITE_EPS ) return m_splines[0].evalNormal(0.0);
    if ( t > m_spline_seg - HERMITE_EPS ) return m_splines[m_spline_seg - 1].evalNormal(1.0);
    return m_splines[static_cast<int>(t)].evalNormal(t - std::floor(t));
}


double HermiteCurve::arcLengthApprox(double t) const
{
    assert(m_spline_seg > 0);
    if ( t < HERMITE_EPS ) return 0.0;
    if ( t > m_spline_seg - HERMITE_EPS ) return m_lens[m_spline_seg - 1];
    int spline_id = static_cast<int>(t);
    return m_splines[spline_id].arcLengthApprox(t - spline_id) + (spline_id ? m_lens[spline_id - 1] : 0.0);
}


double HermiteCurve::arcLengthInvApprox(double len) const
{
    assert(m_spline_seg > 0);
    if ( len < HERMITE_EPS ) return 0.0;
    if ( len > m_lens[m_spline_seg - 1] - HERMITE_EPS ) return static_cast<double>(m_spline_seg);
    int spline_id = findSegId(len);
    if ( spline_id ) len -= m_lens[spline_id - 1];
	return static_cast<double>(spline_id) + m_splines[spline_id].arcLengthInvApprox(len);
}


void HermiteCurve::segLengths(std::vector<double> &length_list) const
{
    assert(m_spline_seg > 0);
    length_list.resize(m_spline_seg);
    length_list[0] = m_lens[0];
	for (int i = 1; i < m_spline_seg; i++)
        length_list[i] = m_lens[i] - m_lens[i - 1];
}


int HermiteCurve::findSegId(double curve_length) const
{
    assert(m_spline_seg > 0);
    assert(curve_length > -HERMITE_EPS && m_lens[m_spline_seg - 1] - curve_length > -HERMITE_EPS);
    if ( curve_length < HERMITE_EPS ) return 0;
    if ( curve_length > m_lens[m_spline_seg - 1] - HERMITE_EPS ) return m_spline_seg - 1;
    return static_cast<int>(std::lower_bound(m_lens.begin(), m_lens.end(), curve_length) - m_lens.begin());
}


void HermiteCurve::output(int n, Eigen::Vector3d *bufferPosition,
    Eigen::Vector3d *bufferTangent, Eigen::Vector3d *bufferNormal) const
{
    assert(m_spline_seg > 0 && n > 1);
    for ( int i = 0; i < n; ++i ) {
        double t = static_cast<double>(i)*m_spline_seg/(n - 1);
        bufferPosition[i] = eval(t);
        if ( bufferTangent ) bufferTangent[i] = evalTangent(t);
        if ( bufferNormal ) bufferNormal[i] = evalNormal(t);
    }
}
