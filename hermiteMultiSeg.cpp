#include <algorithm>
#include "hermiteSeg.h"
#include "hermiteMultiSeg.h"

#define HERMITE_EPS     1e-5
#define POLYSOLVER_EPS  1e-8

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