#pragma once
#ifndef __HERMITE_MULTISEG_SPLINE_H_
#define __HERMITE_MULTISEG_SPLINE_H_

#include <vector>
#include <sstream>
#include <Eigen/Dense>
#include "Util.h"

class HermiteSpline_multiSeg {
public:
	HermiteSpline_multiSeg(const char* filename) {
		init_splines(filename);
	}
	void init_splines(const char* filename);
	Eigen::Vector3d eval(double t);
	Eigen::Vector3d evalTangent(double t);
	Eigen::Vector3d evalCurvature(double t);
	double HermiteSpline_multiSeg::arcLengthInvApprox(double len, int subdiv);

	/* Get number of segments for the multi-seg spline */
	inline int get_seg_num() const {
		return m_spline_seg;
	}

	/* Get seg_i of the multi-seg spline */
	inline HermiteSpline get_spline(int id) {
		return m_splines[id];
	}

	/* Get the list contains length of each segment */
	std::vector<double> segLengths(int subdiv);

	/* Get the total length of the multi-seg spline */
	double totalLength(int subdiv);

	/* Find which segment belongs to this curve-length*/
	int findSegId(double curve_length, int subdiv);

protected:
	std::vector<HermiteSpline> m_splines;
	int m_spline_seg;
};

#endif