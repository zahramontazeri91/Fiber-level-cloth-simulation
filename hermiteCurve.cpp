#include <algorithm>
#include "hermiteCurve.h"


void HermiteCurve::init_norm(const char* pntsFILE, const char* normsFILE, int subdiv) {
    assert(pntsFILE);
	std::ifstream fin(pntsFILE);
    assert(!fin.fail());

    int n;
    fin >> n;

    std::vector<Eigen::Vector3d> pts(n), norms(n);
    for ( int i = 0; i < n; ++i ) fin >> pts[i][0] >> pts[i][1] >> pts[i][2];

    init_norm(pts, subdiv);

	//TODO : to print normals comment buildPlanes and PlanesIntersections2D in CrossSection::init()
	printNormals(normsFILE, subdiv);
}

void HermiteCurve::printNormals(const char* normsFILE, const int subdiv) {
	assert(normsFILE);
	std::ofstream fout(normsFILE);
	fout << m_splines.size() + 1  << '\n';
	for (int i = 0; i < m_splines.size(); ++i) {
		Eigen::Vector3d n = m_splines[i].evalNormal(0);
		fout << n[0] << ' ' << n[1] << ' ' << n[2] << std::endl;
	}
	int i = m_splines.size() - 1;
	Eigen::Vector3d n = m_splines[i].evalNormal(1);
	fout << n[0] << ' ' << n[1] << ' ' << n[2] << std::endl;
	std::cout << "Normals are written to curveNormals.txt \n";
	fout.close();
}

void HermiteCurve::init(const char* pntsFILE, const char* normsFILE, int subdiv) {
	// import the points
	assert(pntsFILE);
	std::ifstream fin(pntsFILE);
	assert(!fin.fail());
	int n;
	fin >> n;
	std::vector<Eigen::Vector3d> pts(n);
	for (int i = 0; i < n; ++i) fin >> pts[i][0] >> pts[i][1] >> pts[i][2];

	// import the normals 
	assert(normsFILE);
	std::ifstream fin2(normsFILE);
	assert(!fin2.fail());
	fin2 >> n;
	std::vector<Eigen::Vector3d> norms(n);
	for (int i = 0; i < n; ++i) fin2 >> norms[i][0] >> norms[i][1] >> norms[i][2];

	init(pts, norms, subdiv);
}


void HermiteCurve::init_norm(const std::vector<Eigen::Vector3d> &pts, int subdiv) //subdiv for each segment
{
    initPoints(pts);
	
    m_splines[0].build(subdiv, m_splines[0].evalPrincipalNormal(0.0));

    for ( int i = 1; i < m_spline_seg; ++i )
        m_splines[i].build(subdiv, m_splines[i - 1].evalNormal(1.0));
	
    m_lens.resize(m_spline_seg);
    for ( int i = 0; i < m_spline_seg; ++i ) {
        m_lens[i] = m_splines[i].totalLength();
        if ( i ) m_lens[i] += m_lens[i - 1];
    }
}


void HermiteCurve::init(const std::vector<Eigen::Vector3d> &pts, const std::vector<Eigen::Vector3d> &norms, int subdiv)
{
    initPoints(pts);
	if (pts.size() != norms.size())
		std::cout << "number of curve points and norms don't match! " << pts.size() << " " << norms.size() << std::endl;
    assert(pts.size() == norms.size());

    for ( int i = 0; i < m_spline_seg; ++i )
        m_splines[i].build(subdiv, norms[i], norms[i + 1]);

    m_lens.resize(m_spline_seg);
    for ( int i = 0; i < m_spline_seg; ++i ) {
        m_lens[i] = m_splines[i].totalLength();
        if ( i ) m_lens[i] += m_lens[i - 1];
    }
}


void HermiteCurve::initPoints(const std::vector<Eigen::Vector3d> &pts)
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

void HermiteCurve::getRotatedFrame(double t, Eigen::Vector3d &ex, Eigen::Vector3d &ey, Eigen::Vector3d &ez) const {
	//rotating Frenet frame by 90 is same as follows:
	Eigen::Vector3d T = evalTangent(t);
	Eigen::Vector3d N = evalNormal(t);
	Eigen::Vector3d B = T.cross(N);
	
	//ez = T;
	//ex = N;
	//ey = B;

	ez = T;
	ex = -1*B;
	ey = N;

	//ex = -1*N;
	//ey = -1*B;

	assert(ex.norm() - 1.f < 1e-5);
	assert(ey.norm() - 1.f < 1e-5);
	assert(ez.norm() - 1.f < 1e-5);
}