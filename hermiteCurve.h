#pragma once
#ifndef __HERMITE_MULTISEG_SPLINE_H_
#define __HERMITE_MULTISEG_SPLINE_H_

#include <vector>
#include <sstream>
#include <Eigen/Dense>
#include "Util.h"
#include "hermiteSeg.h"

class HermiteCurve {
public:
    inline HermiteCurve(const char* filename)               { init_splines(filename); }

	void init_splines(const char* filename);

	Eigen::Vector3d eval(double t) const;
	Eigen::Vector3d evalTangent(double t) const;
	Eigen::Vector3d evalCurvature(double t) const;

    double arcLength(double len) const;
	double arcLengthInvApprox(double len) const;

	/* Get number of segments for the multi-seg spline */
    inline int get_seg_num() const                          { return static_cast<int>(m_splines.size()); }

	/* Get seg_i of the multi-seg spline */
	inline const HermiteSpline &get_spline(int id) const    { return m_splines[id]; }

	/* Get the list contains length of each segment */
	std::vector<double> segLengths();

	/* Get the total length of the multi-seg spline */
	double totalLength();

	/* Find which segment belongs to this curve-length*/
	int findSegId(double curve_length);

protected:
    int m_spline_seg;
	std::vector<HermiteSpline> m_splines;
};

#endif