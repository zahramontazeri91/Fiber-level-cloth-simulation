#pragma once
#include "Util.h"
#include "Curve.h"

typedef std::vector<vec3f> fiber_t;
typedef std::vector<fiber_t> yarn_t;

namespace Fiber {

	/* Define a fiber */
	struct Fiber {
		std::vector<vec3f> vertices;	 // world-space positions for this fiber
		vec3f init_vertex;				 // initial vertex of this fiber on z=0 normal plane, separate from the previous vertices
										 // because we might need fiber migration 
		float init_radius;				 // initial radius sampled with Section 4.1 Cross-Sectional Fiber Distribution
		float init_theta;				 // initial theta for this fiber in circle
		float init_migration_theta;		 // initial migration theta for this fiber
		void clear() {
			vertices.clear();
		}
	};

	/* Define a ply */
	struct Ply {
		/* Fibers in this ply */
		std::vector<Fiber> fibers;

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

								 /* Parameters and functions for Section 4.4 Hairiness */
		float mu, sigma;		 // Flyaway fiber length ~ NormalDistribution(mu, sigma)
		int flyaway_num;		 // Flyaway fiber number in this ply
		float fly_step_size;	 // Flyaway step size 

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

		/* Define a set of parameters needed for compression*/    //TODO: unify this with the other crossSection/Elllipse Struct
		struct compress {
			float z;
			float theta;
			float ellipse_long;
			float ellipse_short;
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

		/* Build and initialize a 1-ply yarn from a given file */
		void build(const char *yarnfile, const int num_ply);

		/* Parse config file from disk */
		void parse(const char *filename = NULL);

		/* Simulate yarn */
		void yarn_simulate();

		/* compress yarn with theta, direction, a and b, cross section ellipse param. */
		void compress_yarn(const char* filename);

		/* Write simulated data (single yarns) to disk */
		void write_yarn(const char* filename);

		/* map the straight yarn on a curved spline */
		void curve_yarn(const char* filename, bool scaleXY = false);
	};
}