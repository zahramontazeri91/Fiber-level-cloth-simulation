#include <algorithm>
#include "hermiteCurve.h"


void HermiteCurve::init_splines(const char* filename) {
	// read the spline parameters file
	std::ifstream fin;
	if (filename != NULL)
		fin.open(filename);

    fin >> m_spline_seg;
    m_splines.resize(m_spline_seg);

	for (int i = 0; i < m_spline_seg; i++) {
        Eigen::Vector3d p0, p1, m0, m1;
        fin >> p0[0] >> p0[1] >> p0[2] >> p1[0] >> p1[1] >> p1[2] >> m0[0] >> m0[1] >> m0[2] >> m1[0] >> m1[1] >> m1[2];

		// initialize the splines
        m_splines[i].init(p0, p1, m0, m1);
	}
}

Eigen::Vector3d HermiteCurve::eval(double t) const {
	int spline_id = int(t); // curve from t = 0 to t = 1 defined by first segment and so on
	if (t == m_spline_seg)
		return m_splines[spline_id - 1].eval(1.0);
	return m_splines[spline_id].eval(t - double(spline_id));
}

Eigen::Vector3d HermiteCurve::evalTangent(double t) const {
	int spline_id = int(t); // curve from t = 0 to t = 1 defined by first segment and so on
	if (t == m_spline_seg)
		return m_splines[spline_id - 1].evalTangent(1.0);
	return m_splines[spline_id].evalTangent(t - double(spline_id));
}

Eigen::Vector3d HermiteCurve::evalCurvature(double t) const {
	int spline_id = int(t);
	if (t == m_spline_seg)
		return m_splines[spline_id - 1].evalCurvature(1.0);
	return m_splines[spline_id].evalCurvature(t - double(spline_id));
}

double HermiteCurve::arcLengthInvApprox(double len) const {

	int spline_id = int(len);
	if (len == m_spline_seg)
		return m_splines[spline_id - 1].arcLengthInvApprox(1.0);
	return m_splines[spline_id].arcLengthInvApprox(len - double(spline_id));
}

std::vector<double> HermiteCurve::segLengths() {
	std::vector<double> length_list;
	for (int i = 0; i < m_spline_seg; i++) {
		double length_seg_i = m_splines[i].totalLength();
		length_list.push_back(length_seg_i);
	}
	return length_list;
}

double HermiteCurve::totalLength() {
	std::vector<double> length_list = segLengths();
	double sum = 0.0;
	for (int i = 0; i < length_list.size(); i++) {
		sum += length_list[i];
	}
	return sum;
}

int HermiteCurve::findSegId(double curve_length) { //TODO:move length_list to class member
	int segId = 0;
	double length_sofar = 0.0;
	std::vector<double> length_list = segLengths();

	for (int i = 0; i < length_list.size(); i++) {
		if ((curve_length - length_sofar) >= 0) {
			length_sofar += length_list[i];
			segId = i;
			continue;
		}
		else {
			segId = i - 1;
			break;
		}
	}
	return segId;
}