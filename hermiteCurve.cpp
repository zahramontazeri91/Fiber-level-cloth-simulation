#include <algorithm>
#include "hermiteCurve.h"


void HermiteCurve::init(const char* pntsFILE, const char* normsFILE, int subdiv) {

    assert(pntsFILE);
	std::ifstream fin(pntsFILE);

	if (fin.fail())
		std::cout << pntsFILE;

    assert(!fin.fail());

    int n;
    fin >> n;

    std::vector<Eigen::Vector3d> pts(n), norms(n);
    for ( int i = 0; i < n; ++i ) fin >> pts[i][0] >> pts[i][1] >> pts[i][2];

	assert(pts.size()>2);

    init(pts, subdiv);
	printNormals(normsFILE, subdiv);
}

// not using this anymore
void HermiteCurve::printBiNormals(const char* binormsFILE, const int subdiv) {
	assert(binormsFILE);
	std::ofstream fout(binormsFILE);
	fout << m_splines.size() + 1 << '\n';
	for (int i = 0; i < m_splines.size(); ++i) {
		Eigen::Vector3d n = m_splines[i].evalNormal(0);
		Eigen::Vector3d t = m_splines[i].evalTangent(0);
		Eigen::Vector3d b = t.cross(n);
		fout << b[0] << ' ' << b[1] << ' ' << b[2] << std::endl;
	}
	int i = m_splines.size() - 1;
	Eigen::Vector3d n = m_splines[i].evalNormal(1);
	Eigen::Vector3d t = m_splines[i].evalTangent(1);
	Eigen::Vector3d b = t.cross(n);
	fout << b[0] << ' ' << b[1] << ' ' << b[2] << std::endl;
	//std::cout << "Normals are written to the file \n";
	fout.close();
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
	//std::cout << "Normals are written to the file \n";
	fout.close();
}

void HermiteCurve::init_norm(const char* pntsFILE, const char* normsFILE, int subdiv) {

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
	for (int i = 0; i < n; ++i) {
		fin2 >> norms[i][0] >> norms[i][1] >> norms[i][2];
		//std::cout << norms[i][0] << norms[i][1] << norms[i][2] << std::endl;
	}

	init_norm(pts, norms, subdiv);
}

void HermiteCurve::init_seg(const char* pntsFILE, const int start, const int end, int subdiv) {

	assert(pntsFILE);
	std::ifstream fin(pntsFILE);
	assert(!fin.fail());

	int n;
	fin >> n;
	assert(start <= end);

	const int seg_size = end - start + 1;
	std::vector<Eigen::Vector3d> allPts(n), pts(seg_size);
	int j = 0;
	for (int i = 0; i < n; ++i) {
		fin >> allPts[i][0] >> allPts[i][1] >> allPts[i][2];
		if (i >= start && i <= end) {
			pts[j] = allPts[i];
			j++;
		}
	}
	assert(pts.size()>2);

	init(pts, subdiv);
	//printNormals(normsFILE, subdiv);
}

void HermiteCurve::init_principleNormal(const char* pntsFILE, const char* binormsFILE, int subdiv) {

	assert(pntsFILE);
	std::ifstream fin(pntsFILE);
	assert(!fin.fail());

	int n;
	fin >> n;

	std::vector<Eigen::Vector3d> pts(n), norms(n);
	for (int i = 0; i < n; ++i) fin >> pts[i][0] >> pts[i][1] >> pts[i][2];
	assert(pts.size()>2);

	init_principleNormal(pts, subdiv);
	printBiNormals(binormsFILE, subdiv);
}

/*mitsuba implementation: given vector a, it return a,b,c such that they form an orthogonal basis*/
void coordinateSystem(const vec3f &a, vec3f &b, vec3f &c) {
	if (std::abs(a.x) > std::abs(a.y)) {
		float invLen = 1.0f / std::sqrt(a.x * a.x + a.z * a.z);
		c = vec3f(a.z * invLen, 0.0f, -a.x * invLen);
	}
	else {
		float invLen = 1.0f / std::sqrt(a.y * a.y + a.z * a.z);
		c = vec3f(0.0f, a.z * invLen, -a.y * invLen);
	}
	b = cross(c, a);
}

#if 0
void HermiteCurve::init(const std::vector<Eigen::Vector3d> &pts, int subdiv) //subdiv for each segment
{
	initPoints(pts);
	m_splines[0].build(subdiv, m_splines[0].evalPrincipalNormal(0.0));
	for (int i = 1; i < m_spline_seg; ++i)
		m_splines[i].build(subdiv, m_splines[i - 1].evalNormal(1.0));

	m_lens.resize(m_spline_seg);
	for (int i = 0; i < m_spline_seg; ++i) {
		m_lens[i] = m_splines[i].totalLength();
		if (i) m_lens[i] += m_lens[i - 1];
	}


}
#endif 

#if 1
void HermiteCurve::init(const std::vector<Eigen::Vector3d> &pts, int subdiv) //subdiv for each segment
{
	initPoints(pts);

	/* Find the first vertex that curvature doesn't vanish */
	int firstIndx = 0;
	bool isFound = false;
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 0; i < m_spline_seg; ++i) {
		Eigen::Vector3d q = m_splines[i].evalCurvature(0.0);
		if (q.norm() > HERMITE_EPS) {
			firstIndx = i;
			isFound = true;
			break;
		}
	}

	//assert(isFound && "Assuming there exist one point in the curve that its curvature doesn't vanish");
	//if (firstIndx>0)
		//std::cout << firstIndx << "-th vertex has non-vanishing curvature. \n";

	if (isFound) {
		m_splines[firstIndx].build(subdiv, m_splines[firstIndx].evalPrincipalNormal(0.0));
	}
	else {//when all the points of the window have vanished curvature
		firstIndx = 0; //pick a fixed normal for the first point
		Eigen::Vector3d t = m_splines[firstIndx].evalTangent(0.0);
		vec3f tang(t(0), t(1), t(2));
		vec3f norm, binorm;
		coordinateSystem(tang, norm, binorm);
		Eigen::Vector3d n;
		n << norm.x, norm.y, norm.z;
		m_splines[firstIndx].build(subdiv, n);

	}

    
	for (int i = firstIndx+1; i < m_spline_seg; ++i) {
		m_splines[i].build(subdiv, m_splines[i - 1].evalNormal(1.0), Eigen::Vector3d::Zero());
	}
	for (int i = firstIndx - 1; i >= 0; --i) { //go backward
		m_splines[i].build(subdiv, Eigen::Vector3d::Zero(), m_splines[i + 1].evalNormal(0.0));
	}
    m_lens.resize(m_spline_seg);
    for ( int i = 0; i < m_spline_seg; ++i ) {
        m_lens[i] = m_splines[i].totalLength();
        if ( i ) m_lens[i] += m_lens[i - 1];
    }

	//DEBUG:
	/*NOTE that tg is not perpendicular to norm for t=0 or t=1 exactly. */
	//Eigen::Vector3d n = m_splines[10].evalNormal(1.0);
	//Eigen::Vector3d t = m_splines[10].evalTangent(1.0);
	//std::cout << "hermitecurve build " << n.dot(t) << std::endl;

	//Eigen::Vector3d t0 = m_splines[10].evalTangent(1.0);
	//Eigen::Vector3d t1 = m_splines[11].evalTangent(0.0);
	//std::cout << "hermitecurve build \n " << t0 << std::endl << t1 << std::endl;
}
#endif

void HermiteCurve::init_norm(const std::vector<Eigen::Vector3d> &pts, const std::vector<Eigen::Vector3d> &norms, int subdiv)
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

void HermiteCurve::init_principleNormal(const std::vector<Eigen::Vector3d> &pts, int subdiv) //subdiv for each segment
{
	initPoints(pts);
	
	//Eigen::Vector3d q = m_splines[0].evalCurvature(0.0);
	//Eigen::Vector3d v = m_splines[0].evalTangent(0.0);
	//std::cout << q.norm() << " " << v.norm() << " " << v.cross(q).cross(v) << " \n (init_principleNormal) \n";
	//assert(q.norm() > HERMITE_EPS);
	//assert(v.norm() > HERMITE_EPS);
	//assert(v.cross(q).cross(v).norm() > HERMITE_EPS);


	m_splines[0].build(subdiv, m_splines[0].evalPrincipalNormal(0.0));
	for (int i = 1; i < m_spline_seg; ++i) {
		// norm0: should be t=1 defines the very last segment of the spline[i-1] and for norm1: t = 1 as it will be interpolated in build()
		m_splines[i].build(subdiv, m_splines[i - 1].evalPrincipalNormal(1.0), m_splines[i].evalPrincipalNormal(1.0));
	}

	m_lens.resize(m_spline_seg);
	for (int i = 0; i < m_spline_seg; ++i) {
		m_lens[i] = m_splines[i].totalLength();
		if (i) m_lens[i] += m_lens[i - 1];
	}
}

void HermiteCurve::initPoints(const std::vector<Eigen::Vector3d> &pts)
{
    assert(pts.size()>2);

    m_spline_seg = static_cast<int>(pts.size()) - 1;
    m_splines.resize(m_spline_seg);
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
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

void HermiteCurve::assign_upsample(std::vector<Eigen::Vector3d> &all_pts, std::vector<Eigen::Vector3d> &all_tg, std::vector<Eigen::Vector3d> &all_norm) {
	all_pts.clear();
	all_tg.clear();
	all_norm.clear();

	/* upsample a curve by 2 */
	float t0 = 0.0;
	float t1 = 0.5;

	//duplicate the last point 
	int i = 0;
	Eigen::Vector3d pnt = m_splines[i].eval(0.0);
	Eigen::Vector3d tg = m_splines[i].evalTangent(0.0);
	Eigen::Vector3d norm = m_splines[i].evalNormal(0.0);
	all_pts.push_back(pnt);
	all_tg.push_back(tg);
	all_norm.push_back(norm);

	pnt = m_splines[i].eval(0.3);
	tg = m_splines[i].evalTangent(0.3);
	norm = m_splines[i].evalNormal(0.3);
	///assert(std::abs(norm.dot(tg)) < eps); ## TODO
	all_pts.push_back(pnt);
	all_tg.push_back(tg);
	all_norm.push_back(norm);

	pnt = m_splines[i].eval(0.6);
	tg = m_splines[i].evalTangent(0.6);
	norm = m_splines[i].evalNormal(0.6);
	///assert(std::abs(norm.dot(tg)) < eps); ## TODO
	all_pts.push_back(pnt);
	all_tg.push_back(tg);
	all_norm.push_back(norm);

	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 1; i < m_splines.size(); ++i) {
		pnt = m_splines[i].eval(t0);
		tg = m_splines[i].evalTangent(t0);
		norm = m_splines[i].evalNormal(t0);
		///assert(std::abs(norm.dot(tg)) < eps); ## TODO
		all_pts.push_back(pnt);
		all_tg.push_back(tg);
		all_norm.push_back(norm);

		pnt = m_splines[i].eval(t1);
		tg = m_splines[i].evalTangent(t1);
		norm = m_splines[i].evalNormal(t1);
		///assert(std::abs(norm.dot(tg)) < eps); ## TODO
		all_pts.push_back(pnt);
		all_tg.push_back(tg);
		all_norm.push_back(norm);

	}
	i = m_splines.size() - 1;
	pnt = m_splines[i].eval(1);
	tg = m_splines[i].evalTangent(1);
	norm = m_splines[i].evalNormal(1);
	all_pts.push_back(pnt);
	all_tg.push_back(tg);
	all_norm.push_back(norm);

}

void HermiteCurve::assign(std::vector<Eigen::Vector3d> &all_pts, std::vector<Eigen::Vector3d> &all_tg, std::vector<Eigen::Vector3d> &all_norm) {
	all_pts.clear();
	all_tg.clear();
	all_norm.clear();

	float t = 0.0;

	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 0; i < m_splines.size(); ++i) {
		Eigen::Vector3d pnt = m_splines[i].eval(t);
		Eigen::Vector3d tg = m_splines[i].evalTangent(t);
		Eigen::Vector3d norm = m_splines[i].evalNormal(t);
		//std::cout << norm.dot(tg) << std::endl;
		///assert(std::abs(norm.dot(tg)) < eps); ## TODO
		all_pts.push_back(pnt);
		all_tg.push_back(tg);
		all_norm.push_back(norm);
	}
	int i = m_splines.size() - 1;
	Eigen::Vector3d pnt = m_splines[i].eval(1);
	Eigen::Vector3d tg = m_splines[i].evalTangent(1);
	Eigen::Vector3d norm = m_splines[i].evalNormal(1);
	all_pts.push_back(pnt);
	all_tg.push_back(tg);
	all_norm.push_back(norm);
}

void HermiteCurve::assign_twist(const char* twistFile, std::vector<Eigen::Vector3d> &all_pts, std::vector<Eigen::Vector3d> &all_tg, 
	std::vector<Eigen::Vector3d> &all_norm_rot, const int upsample) {
	all_pts.clear();
	all_tg.clear();
	all_norm_rot.clear();
	std::vector<Eigen::Vector3d> all_norm;

	float t0 = 0.0;
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 0; i < m_splines.size(); ++i) {
		Eigen::Vector3d pnt = m_splines[i].eval(t0);
		Eigen::Vector3d tg = m_splines[i].evalTangent(t0);
		Eigen::Vector3d norm = m_splines[i].evalNormal(t0);
		all_pts.push_back(pnt);
		all_tg.push_back(tg);
		all_norm.push_back(norm);
		if (i == 0) {
			for (int d = 1; d < 2* upsample - 1; d++) { // for upsampling, one sample will already counted for t=0
				float t = float(d) * (1.0 / (2.0*upsample - 1) );
				pnt = m_splines[i].eval(t);
				tg = m_splines[i].evalTangent(t);
				norm = m_splines[i].evalNormal(t);
				all_pts.push_back(pnt);
				all_tg.push_back(tg);
				all_norm.push_back(norm);
			}
			continue;
		}
		for (int d = 1; d < upsample; d++) {
			float t = float(d) * ( 1.0 / upsample);
			pnt = m_splines[i].eval(t);
			tg = m_splines[i].evalTangent(t);
			norm = m_splines[i].evalNormal(t);
			all_pts.push_back(pnt);
			all_tg.push_back(tg);
			all_norm.push_back(norm);
		}
	}
	int i = m_splines.size() - 1 ; //last segment
	Eigen::Vector3d pnt = m_splines[i].eval(1);
	Eigen::Vector3d tg = m_splines[i].evalTangent(1);
	Eigen::Vector3d norm = m_splines[i].evalNormal(1);
	all_pts.push_back(pnt);
	all_tg.push_back(tg);
	all_norm.push_back(norm);

	/* twist the normals by the twisting angles */
	//all_norm_rot = all_norm;
	twistNormals(twistFile, all_tg, all_norm, all_norm_rot);
}
void HermiteCurve::assign_twist(const std::vector<float> &twists, std::vector<Eigen::Vector3d> &all_pts, std::vector<Eigen::Vector3d> &all_tg,
	std::vector<Eigen::Vector3d> &all_norm_rot, const int upsample) {
	all_pts.clear();
	all_tg.clear();
	all_norm_rot.clear();
	std::vector<Eigen::Vector3d> all_norm;

	float t0 = 0.0;
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 0; i < m_splines.size(); ++i) {
		Eigen::Vector3d pnt = m_splines[i].eval(t0);
		Eigen::Vector3d tg = m_splines[i].evalTangent(t0);
		Eigen::Vector3d norm = m_splines[i].evalNormal(t0);
		all_pts.push_back(pnt);
		all_tg.push_back(tg);
		all_norm.push_back(norm);
		if (i == 0) {
			for (int d = 1; d < 2 * upsample - 1; d++) { // for upsampling, one sample will already counted for t=0
				float t = float(d) * (1.0 / (2.0*upsample - 1));
				pnt = m_splines[i].eval(t);
				tg = m_splines[i].evalTangent(t);
				norm = m_splines[i].evalNormal(t);
				all_pts.push_back(pnt);
				all_tg.push_back(tg);
				all_norm.push_back(norm);
			}
			continue;
		}
		for (int d = 1; d < upsample; d++) {
			float t = float(d) * (1.0 / upsample);
			pnt = m_splines[i].eval(t);
			tg = m_splines[i].evalTangent(t);
			norm = m_splines[i].evalNormal(t);
			all_pts.push_back(pnt);
			all_tg.push_back(tg);
			all_norm.push_back(norm);
		}
	}
	int i = m_splines.size() - 1; //last segment
	Eigen::Vector3d pnt = m_splines[i].eval(1);
	Eigen::Vector3d tg = m_splines[i].evalTangent(1);
	Eigen::Vector3d norm = m_splines[i].evalNormal(1);
	all_pts.push_back(pnt);
	all_tg.push_back(tg);
	all_norm.push_back(norm);

	/* twist the normals by the angles */
	twistNormals(twists, all_tg, all_norm, all_norm_rot);
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
	Eigen::Vector3d n = m_splines[static_cast<int>(t)].evalNormal(t - std::floor(t));

	//assert(std::abs(n.norm() - 1.0) < HERMITE_EPS);
    return n;
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

	ez = T;
	//ex = N;
	//ey = B;
	ex = B;
	ey = N;

	if (ex.norm() - 1.f > 1e-5 || ey.norm() - 1.f > 1e-5 || ez.norm() - 1.f > 1e-5)
		std::cout << ex.norm() - 1.f << " " << ey.norm() - 1.f << " " << ez.norm() - 1.f << std::endl;
	assert(ex.norm() - 1.f < 1e-5);
	assert(ey.norm() - 1.f < 1e-5);
	assert(ez.norm() - 1.f < 1e-5);
}

Eigen::Vector3d HermiteCurve::rotVec3(const float angle, const Eigen::Vector3d &axis, const Eigen::Vector3d &vec) {
	//using "rotation about an arbitrary axis"
	const float C = cos(angle);
	const float S = sin(angle);
	const float t = 1.0 - cos(angle);
	const float ux = axis[0];
	const float uy = axis[1];
	const float uz = axis[2];

	Eigen::Matrix3d M;
	M << t*ux*ux + C, t*ux*uy - S*uz, t*ux*uz + S*uy,
		t*ux*uy + S*uz, t*uy*uy + C, t*uy*uz - S*ux,
		t*ux*uz - S*uy, t*uy*uz + S*ux, t*uz*uz + C;

	return  M*vec;
}

void HermiteCurve::twistNormals(const char* twistFile, const std::vector<Eigen::Vector3d> &all_tg, const std::vector<Eigen::Vector3d> &all_norm, std::vector<Eigen::Vector3d> &all_norm_rot) {

	std::ifstream fin(twistFile);
	if (!fin.is_open())
		std::cout << twistFile << std::endl;
	assert(fin.is_open() && "twist-file is not found!\n");
	int n=0;
	fin >> n;

	//std::cout << "twistNormals " << all_tg.size() << " " << all_norm.size() << std::endl;
	assert(n == all_tg.size() && n == all_norm.size());
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 0; i < n; i++) {

		float twist;
		fin >> twist;

		const Eigen::Vector3d  norm_rot = rotVec3(twist, all_tg[i], all_norm[i]);
		all_norm_rot.push_back(norm_rot);

	}
}

void HermiteCurve::twistNormals(const std::vector<float> &twists, const std::vector<Eigen::Vector3d> &all_tg,
	const std::vector<Eigen::Vector3d> &all_norm, std::vector<Eigen::Vector3d> &all_norm_rot) {

	int n = twists.size();
	if (n != all_tg.size() || n != all_norm.size())
		std::cout << "n: " << n << " all_tg: " << all_tg.size() << " all_norm: " << all_norm.size();
	assert(n == all_tg.size() && n == all_norm.size());
	const int num_of_cores = omp_get_num_procs();
#pragma omp parallel for num_threads(num_of_cores)
	for (int i = 0; i < n; i++) {

		float twist;
		twist = twists[i];

		const Eigen::Vector3d  norm_rot = rotVec3(twist, all_tg[i], all_norm[i]);
		all_norm_rot.push_back(norm_rot);
	}
}