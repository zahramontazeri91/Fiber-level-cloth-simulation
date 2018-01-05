#pragma once
#include "Util.h"
#include "Curve.h"
#include <Eigen/Dense>

//#define IMPROVED_FLYAWAYS


typedef std::vector<vec3f> fiber_t;
typedef std::vector<fiber_t> yarn_t;

typedef std::vector<vec3f> plyIntersect;				     //Plane intersection with each of the plys
typedef std::vector<plyIntersect> yarnIntersect;		     //Plane intersection with whole yarn
typedef std::vector<vec2f> plyIntersect2D;				     //Plane intersection with each of the plys in 2D
typedef std::vector<plyIntersect2D> yarnIntersect2D;		     //Plane intersection with whole yarn in 2D

namespace Fiber {

	/* Define a fiber */
	struct Fiber {
		std::vector<vec3f> vertices;	 // world-space positions for this fiber
#ifndef IMPROVED_FLYAWAYS
		std::vector<std::vector<vec3f> > fly_vertices_list; // world-space positions for fly away part
#endif
		vec3f init_vertex;				 // initial vertex of this fiber on z=0 normal plane, separate from the previous vertices
										 // because we might need fiber migration 
		float init_radius;				 // initial radius sampled with Section 4.1 Cross-Sectional Fiber Distribution
		float init_theta;				 // initial theta for this fiber in circle
		float init_migration_theta;		 // initial migration theta for this fiber
		void clear() {
			vertices.clear();
#ifndef IMPROVED_FLYAWAYS
			for (int i = 0; i < fly_vertices_list.size(); i++)
				fly_vertices_list[i].clear();
			fly_vertices_list.clear();
#endif
		}
	};

	/* Define a ply */
	struct Ply {
		/* Fibers in this ply */
		std::vector<Fiber> fibers;
#ifndef IMPROVED_FLYAWAYS
		int fly_fiber_num;
#endif

		/* Base center and theta on z=0 normal plane, and default radius = yarn_radius/2 */
		vec3f base_center;
		float base_theta;

		/* Parameters and functions for Section 4.1 Cross-Sectional Fiber Distribution */
		float epsilon;
		float R_max;
		float beta;
		float fiberDistribution(float R);
		float sampleR();

		/* Parameters and functions for Section 4.2 Twist/Fiber Migration */
		float alpha;
		float s_i;
		float rho_min, rho_max;
		float helixRadius(float init_r, float init_migration_theta, float theta, bool use_migration);
		void helixXYZ(float init_r, float init_theta, float theta, bool use_migration, float init_migration_theta, float &x, float &y/*, float &z*/);

		/* Parameters and functions for Section 4.3 Two-Ply Yarns/Strand Compression */
		float ellipse_long, ellipse_short; // e_i and d_i in original paper
		bool clock_wise;		 // This parameter determines twisting direction of the fibers in ply 

#ifdef IMPROVED_FLYAWAYS
								 /* Parameters for Fly-aways */
		float flyaway_hair_density;
		float flyaway_hair_ze_mu, flyaway_hair_ze_sigma;
		float flyaway_hair_r0_mu, flyaway_hair_r0_sigma;
		float flyaway_hair_re_mu, flyaway_hair_re_sigma;
		float flyaway_hair_pe_mu, flyaway_hair_pe_sigma;

		float flyaway_loop_density;
		float flyaway_loop_r1_mu, flyaway_loop_r1_sigma;
#else
								 /* Parameters and functions for Section 4.4 Hairiness */
		float mu, sigma;		 // Flyaway fiber length ~ NormalDistribution(mu, sigma)
		int flyaway_num;		 // Flyaway fiber number in this ply
		float fly_step_size;	 // Flyaway step size 
#endif

	};

	class Yarn {
		/* Plys in this yarn */
		//std::vector<Ply> plys; make this public

		bool use_migration;		// Fiber migration
		bool clock_wise;		// This parameter determines twisting direction of the plys in yarn 
		float z_step_size;		// Yarn grows along z-axis, this parameter determines the growth step size
		float z_step_num;		// Yarn grows with given step number

								/* Coaxial helix model: x = R(theta)cos(theta), y = R(theta)sin(theta), z = alpha * theta / (2*PI) */
		float yarn_alpha;		// Plys twisting in yarn
		float yarn_radius;		// Radius for yarn circle in which plys twist

#ifdef IMPROVED_FLYAWAYS
		bool use_flyaways;
#endif

		std::string config_file;  // Config parameters filename
		std::string output_file;  // Data output filename 

		Curve z_curve;			  // Fit yarn along the given curve 

		AABB aabb_procedural;	  // Bounding box for all procedural vertices (no flyaway part)
		AABB aabb_micro_ct;		  // Bounding box for Micro-CT volume data

		omp_lock_t lock;		  // Perlin noise generator needs it

		float compress_theta;     // Angle that the cross section should be rotated to apply the yarn compression (radian)
		float compress_long, compress_short;

	public:
		std::vector<Ply> plys;

		struct Transform {
			Eigen::Matrix2f R;
			Eigen::Matrix2f S;
		};

		/* Define a set of parameters needed for compression*/    //TODO: unify this with the other crossSection/Elllipse Struct
		struct Compress {
			float ellipse_long;
			float ellipse_short;
			float ellipse_theta;
			float rotation;
		};

		struct CenterLine {
			//model each segment of yarn with one cycle sinisoid a*sin(x+b) + c
			float a;
			float b;
			float c;
			float d;
		};

		void setStepSize(const float ss) {
			this->z_step_size = ss;
		}
		void setStepNum(const int sn) {
			this->z_step_num = sn;
		}
		float getStepSize() const {
			return this->z_step_size;
		}
		int getStepNum() const {
			return this->z_step_num;
		}
		int getPlyNum() const {
			return this->plys.size();
		}
		float getYarnRadius() const {
			return this->yarn_radius;
		}
		float scaleFacFromRadius(
			const float newYR, const int newN, const float newSS, const float meshScaler
		)
		{
			std::cout << "scaleFacFromRadius()" << std::endl;
			float oldYR = this->yarn_radius;
			float oldSS = this->z_step_size;
			int oldN = this->z_step_num;
			float oldLen = oldSS * oldN;

			float newLen = newSS * newN;

			float scaleLen = newLen / oldLen;

			std::cout << "oldLen = " << oldLen << " newLen = " << newLen << " scaleLen=" << scaleLen << std::endl;


			this->aabb_micro_ct.pMin.z *= scaleLen;
			this->aabb_micro_ct.pMax.z *= scaleLen;

			const float THE_NEW_SS = 13;

			this->setStepSize(THE_NEW_SS);
			this->setStepNum(newLen / THE_NEW_SS);

			std::cout << "SS=" << this->getStepSize() << " SN=" << this->getStepNum() << std::endl;

			this->yarn_alpha *= 4; //   twisting
			this->yarn_alpha *= meshScaler;
			for (int i = 0; i < plys.size(); i++) {
				this->plys[i].alpha *= 4; //  twisting
				this->plys[i].alpha *= meshScaler;
			}

			const float scalingFactor = newYR / oldYR;
			return scalingFactor;

		}

	public:
		Yarn();
		~Yarn();

		/* Build and initialize a yarn from a given file */
		void build(const char *yarnfile, const int num_ply);
		void yarnCenter(const char *yarnfile, const char *yarnCenterfile);

		/* Parse config file from disk */
		void parse(const char *filename = NULL);

		/* Simulate yarn */
		void assignParameterizePlyCenters(const char *plyCenterFile);
		void assignPlyCenters(const char *plyCenterFile);
		float extractFiberTwist(const char *fiberTwistFile);
		void yarn_simulate(const char *plyCenter, const char *fiberTwistFile);
		void yarn_simulate();

		/* Simulate yarn with flyaways*/
		/* Simulate ply  */
		void simulate_ply();
		/* Load unrolled plys*/
		void roll_plys(const int K, const std::string &ply_fn, const std::string &fiber_fn);
		/* Write simulated data (separate plys) to disk */
		void write_plys(const char *filename = NULL);

		/* compress yarn with theta, direction, a and b, cross section ellipse param. */
		void readCompressFile(const char* compress_R, const char* compress_S, std::vector<Transform> &all_Transform);
		/* Find the spanning circle for generated yarn before applying the compression */
		void fitCircle(const yarnIntersect2D &pts, float &radius);
		void yarn2crossSections(std::vector<yarnIntersect2D> &itsLists);
		void getPlyCenter(std::vector<std::vector<vec2f>> &plycenters);
		void compress_yarn(const char* compress_R, const char* compress_S);
		void compress_yarn(std::vector<std::vector<Eigen::MatrixXf>> &all_mat_S, std::vector<std::vector<float>> &all_theta_R,
			std::vector<std::vector<Eigen::MatrixXf>> &all_T);
		void compress_yarn3D(const char* deformGrad);

		/* Write simulated data (single yarns) to disk */
		void write_yarn(const char* filename);

		/* map the straight yarn on a curved spline */
		void curve_yarn(const char* pntsFile, const char* normsFile = "", bool scaleXY = false);

		void plotIntersections(const char* filename, const float trimPercent);

		/* Find the yarn-radius for each cross-section */
		//void shapeCrossSection(yarnIntersect2D &its, float &rLong, float &rShort);
	};
}