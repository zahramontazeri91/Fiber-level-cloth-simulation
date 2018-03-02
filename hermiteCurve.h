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
    inline HermiteCurve() {
        m_spline_seg = 0;
    }

	// init1: 
	void init(const char* pntsFILE, const char* normsFILE, int subdiv = 10);
    void init(const std::vector<Eigen::Vector3d> &pts, int subdiv = 10);

	//initi2: given normfile
	void init_norm(const char* pntsFILE, const char* normsFILE, int subdiv = 10);
    void init_norm(const std::vector<Eigen::Vector3d> &pts, const std::vector<Eigen::Vector3d> &norms, int subdiv = 10);

	//initi3: generate principle normal for the whole curve
	void init_principleNormal(const char* pntsFILE, const char* normsFILE, int subdiv);
	void init_principleNormal(const std::vector<Eigen::Vector3d> &pts, int subdiv);

	//init4: generate normals for a segment of the curve [a,b]
	void init_seg(const char* pntsFILE, const int start, const int end, int subdiv);


	Eigen::Vector3d eval(double t) const;
	Eigen::Vector3d evalTangent(double t, bool normalize = true) const;
	Eigen::Vector3d evalCurvature(double t) const;
    Eigen::Vector3d evalNormal(double t) const;

    double arcLengthApprox(double t) const;
	double arcLengthInvApprox(double len) const;

    void output(int n, Eigen::Vector3d *bufferPosition,
        Eigen::Vector3d *bufferTangent = NULL, Eigen::Vector3d *bufferNormal = NULL) const;

	/* Get number of segments for the multi-seg spline */
    inline int get_seg_num() const {
        return m_spline_seg;
    }

	/* Get seg_i of the multi-seg spline */
	inline const HermiteSpline &get_spline(int id) const {
        assert(id >= 0 && id < m_spline_seg);
        return m_splines[id];
    }

	/* Get the list contains length of each segment */
	void segLengths(std::vector<double> &length_list) const;

	/* Get the total length of the multi-seg spline */
    inline double totalLength() const {
        assert(m_spline_seg > 0);
        return m_lens[m_spline_seg - 1];
    }

	/* Find which segment belongs to this curve-length*/
	int findSegId(double curve_length) const;

	/* Returns reference Frenet frame by -90 (so BNT to NBT) */
	void getRotatedFrame(double t, Eigen::Vector3d &ex, Eigen::Vector3d &ey, Eigen::Vector3d &ez) const;

	/* print normals to a file */
	void HermiteCurve::printNormals(const char* normsFILE, const int subdiv = 10);
	void HermiteCurve::printBiNormals(const char* binormsFILE, const int subdiv = 10);

	/* assing points, tangents and normals for specific number of points on the curve */
	void assign(std::vector<Eigen::Vector3d> &all_pts, std::vector<Eigen::Vector3d> &all_tg, std::vector<Eigen::Vector3d> &all_norm);
	void assign_upsample(std::vector<Eigen::Vector3d> &all_pts, std::vector<Eigen::Vector3d> &all_tg, std::vector<Eigen::Vector3d> &all_norm);
protected:
    void initPoints(const std::vector<Eigen::Vector3d> &pts);

    int m_spline_seg;
	std::vector<HermiteSpline> m_splines;

    std::vector<double> m_lens;
};

#endif