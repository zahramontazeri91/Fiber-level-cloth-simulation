#pragma once
#include "nvVector.h"
#include "nvMatrix.h"
#include "Rng.h"
//#include "PerlinNoise.h"
#include <random>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <limits>
#include <time.h>
#include <iomanip>  
#include <omp.h>

const float e = 2.71828182845904523536f;
const float pi = 3.14159265358979323846f;
const float eps = 1e-5f;
const unsigned long long giga = 1024ULL*1024ULL*1024ULL;

typedef nv::vec2<int>		   vec2i;
typedef nv::vec2<float>		   vec2f;
typedef nv::vec3<float>		   vec3f;
typedef nv::vec3<int>          vec3i;
typedef nv::vec3<unsigned int> vec3ui;
typedef nv::vec4<float>        vec4f;
typedef nv::matrix4<float>     mat4f;

#define BINARY_FILE
extern std::string WORK_PATH; 
extern std::vector<Rng> rngs;

static void wait(int seconds) {
	int endwait;
	endwait = clock() + seconds * CLOCKS_PER_SEC;
	while (clock() < endwait) {
		// NO-OP
	}
}

static float rand01() {
	int threadID = omp_get_thread_num();
	assert(threadID >= 0 && threadID < rngs.size());
	return rngs[threadID].rand(0.f, 1.f);
}

static float rand_range(float r_min, float r_max) {
	int threadID = omp_get_thread_num();
	assert(threadID >= 0 && threadID < rngs.size());
	return rngs[threadID].rand(r_min, r_max);
}

static float rand_range(const vec2f &range) {
	return rand_range(range.x, range.y);
}

static float my_round(float number)
{
	return number < 0.0f ? std::ceil(number - 0.5f) : std::floor(number + 0.5f);
}

static bool coin_flip() {
	return rand01() > .5f;
}

/* rotate the vector by angle in radian*/
static vec2f rot2D(const vec2f &v, const float angle) {
	vec2f rotated;
	rotated.x = v.x * std::cos(angle) - v.y * std::sin(angle);
	rotated.y = v.x * std::sin(angle) + v.y * std::cos(angle);
	return rotated;
}

static void print_vec3f(const vec3f &v, const char *name = NULL) {
	if (name) {
		std::cout << name << " " << v.x << " " << v.y << " " << v.z << std::endl;
	} else
		std::cout << v.x << " " << v.y << " " << v.z << std::endl;
}

static void print_vec3f(const vec3f &v, int precision, const char *name = NULL) {
	if (name) {
		std::cout << std::setprecision(precision) << name << " " << v.x << " " << v.y << " " << v.z << std::endl;
	} else
		std::cout << std::setprecision(precision) << v.x << " " << v.y << " " << v.z << std::endl;
}

static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}


static std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

static float clamp(const float &value, const float &low, const float &high) {
	return value < low ? low : (value > high ? high : value); 
}

static struct AABB {
	vec3f pMin;
	vec3f pMax;
	AABB() {
		float float_min = std::numeric_limits<float>::min(),
			  float_max = std::numeric_limits<float>::max();
		pMin = vec3f(float_max, float_max, float_max);
		pMax = vec3f(float_min, float_min, float_min);
	}

	void reset() {
		float float_min = std::numeric_limits<float>::min(),
			  float_max = std::numeric_limits<float>::max();
		pMin = vec3f(float_max, float_max, float_max);
		pMax = vec3f(float_min, float_min, float_min);
	}

	void grow(const vec3f &p) {
		for (int i = 0; i < 3; i++) {
			pMin[i] = std::min(pMin[i], p[i]);
			pMax[i] = std::max(pMax[i], p[i]);
		}
	}
	
	void scale(const float scaler) {
		pMin *= scaler;
		pMax *= scaler;
	}

	bool in(const vec3f &p) const {
		for (int i = 0; i < 3; i++)
			if (p[i] < pMin[i] || p[i] > pMax[i])
				return false;
		return true;
	}

	bool out(const vec3f &p) const {
		return !in(p);
	}

	vec3f shift(const vec3f &p, const float d) const {
		/* shift out vertex with distance d*/
		vec3f p_min = vec3f(vec2f(pMin), 0), p_max = vec3f(vec2f(pMax), 0);
		float r = std::max(nv::length(p_min), nv::length(p_max));
		vec3f p2 = vec3f(vec2f(nv::normalize(vec3f(vec2f(p), 0)) * (r+d)), p.z);
		return p2;
	}

	vec3f lerp(const vec3f &p) const {
		vec3f lerp_p;
		for (int i = 0; i < 3; i++) {
			lerp_p[i] = (p[i] - pMin[i]) / (pMax[i] - pMin[i]); //clamp((p[i] - pMin[i]) / (pMax[i] - pMin[i]), 0.f, 1.f);
		}
		return lerp_p;
	}
};

// Matrix rotation
static vec4f rotvec(const vec3f &a, const vec3f &b) {
	vec3f an = normalize(a);
	vec3f bn = normalize(b);
	vec3f ax = normalize(cross(an, bn));
	float angle = std::acos(std::min(dot(an, bn), 1.f));

	if (length(ax) < eps) {
		vec3f absa = abs(an);
		int index = absa[0] < absa[1] ? 
			(absa[0] < absa[2] ? 0 : 2) :
			(absa[1] < absa[2] ? 1 : 2);
		vec3f c = vec3f(0.f, 0.f, 0.f);
		c[index] = 1.f;
		ax = normalize(cross(an, c));
	}
	return vec4f(ax, angle);
}

static mat4f rotvec2mat(const vec4f &r) {
	const vec3f &ax = vec3f(r);
	const float angle = r.w;
	float s = std::sin(angle);
	float c = std::cos(angle);
	float t = 1 - c;

	vec3f n = normalize(ax);
	float x = n.x, y = n.y, z = n.z;
	mat4f m(
		t*x*x + c,    t*x*y - s*z,  t*x*z + s*y, 0.f,
		t*x*y + s*z,  t*y*y + c,    t*y*z - s*x, 0.f,
		t*x*z - s*y,  t*y*z + s*x,  t*z*z + c,   0.f,
		0.f,		  0.f,			0.f,		 1.f
	);
	return m;
}

static vec3f apply_rotmat(const mat4f &m, const vec3f &a) {
	vec3f b = vec3f(
		dot(vec3f(m.get_column(0)), a),
		dot(vec3f(m.get_column(1)), a),
		dot(vec3f(m.get_column(2)), a)
	);
	return b;
}

static bool BAD(const float &v) {
	return _isnan(v) || !_finite(v);	
}

template <typename T>
static std::vector<size_t> sort_indexes(const std::vector<T> &v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	return idx;
}

template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
	const std::vector<T>& vec,
	Compare& compare)
{
	std::vector<std::size_t> p(vec.size());
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(),
		[&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
	return p;
}

template <typename T>
std::vector<T> apply_permutation(
	const std::vector<T>& vec,
	const std::vector<std::size_t>& p)
{
	std::vector<T> sorted_vec(p.size());
	std::transform(p.begin(), p.end(), sorted_vec.begin(),
		[&](std::size_t i){ return vec[i]; });
	return sorted_vec;
}

template <typename T>
static T statistic_avg(const std::vector<T> &v, int n) {
	const int N = v.size();
	T sum = 0;
	for (size_t i = 0; i != N; ++i)
		sum += v[i];
	T mean = sum / N;
	
	T sqSum = 0;
	for (size_t i = 0; i != N; ++i)
		sqSum += (v[i] - mean) * (v[i] - mean);
	T sigma = std::sqrt(sqSum / (N-1));
	
	T stat_sum = 0;    int count = 0;
	for (size_t i = 0; i != N; ++i)
		if (v[i] < mean + n * sigma && v[i] > mean - n * sigma) {
			stat_sum += v[i];
			count++;
		}
	T stat_avg = stat_sum / count;
	return stat_avg;
}