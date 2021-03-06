#include "Fiber.h"
#include "hermiteCurve.h"
#include "crossSection.h"
#include "fitting.h"

#if 0
#define PERTURB_FIBERS
#ifdef PERTURB_FIBERS
#   define PERTURB_FIBER_PROB          0.9f
#   define PERTURB_FIBER_RATIO         0.25f
#   define PERTURB_FIBER_SMOOTHING     3
#endif
#endif

namespace Fiber {

	float Ply::fiberDistribution(float R) {
		float eTerm = (e - std::powf(e, R / R_max)) / (e - 1);
		float pR = (1 - 2 * epsilon) * std::powf(eTerm, beta) + epsilon;
		return pR;
	}
	float Ply::sampleR() {
		while (true) {
			/* [Attention]: Uniformly sample in circle, r = sqrt(rnd1), theta = 2 * pi * rnd2 */
			/*              This way dS = rdr * dtheta = rnd1^0.5 * 0.5 * rnd1^(-0.5) drnd1 * 2 * pi drnd2 */
			/*						    = pi * drnd1 * drnd2, which is uniform distributed. */
			float radius = std::sqrt((float)rand() / (RAND_MAX)) * R_max, pdf = (float)rand() / (RAND_MAX);
			if (pdf < fiberDistribution(radius))
				return radius;
		}
	}
	float Ply::helixRadius(float init_r, float init_migration_theta, float theta, bool use_migration) { //equation 3 of paper
		float r = init_r;
		if (use_migration)
			r = rho_min * init_r + (rho_max * init_r - rho_min * init_r) * 0.5f * (std::cosf(s_i * (theta)+init_migration_theta) + 1);
		return r;
	}
	void Ply::helixXYZ(float init_r, float init_theta, float theta, bool use_migration, float init_migration_theta, float &x, float &y/*, float &z*/) {
		float r = helixRadius(init_r, init_migration_theta, theta, use_migration);
		x = r * std::cosf(theta + init_theta),
			y = r * std::sinf(theta + init_theta);
		/*z = alpha / (2 * pi) * theta;*/
	}

	void Yarn::build(const char *yarnfile, const int ply_num) {
		//all other parameters already initialized after parse()
		std::ifstream fin (yarnfile);
		if (!fin.is_open()) std::cout << yarnfile << std::endl;
		assert(fin.is_open() && "yarnfile doesn't exist!\n");

		const int num_of_cores = omp_get_num_procs();
		int fiber_num = 0;
		int vrtx_num = 0;
		std::string line;
		std::getline(fin, line);
		fiber_num = atoi(line.c_str()) / ply_num;
		this->plys.resize(ply_num);

		for (int p = 0; p < ply_num; ++p) {
			this->plys[p].fibers.resize(fiber_num);
			for (int f = 0; f < fiber_num; ++f) {
				Fiber &fiber = this->plys[p].fibers[f];
				fiber.clear(); //clear the vertices list 
				std::getline(fin, line);
				vrtx_num = atoi(line.c_str());

				for (int v = 0; v < vrtx_num; ++v) {
					std::getline(fin, line);
					std::vector<std::string> splits = split(line, ' ');
					vec3f vrtx(atof(splits[0].c_str()), atof(splits[1].c_str()), atof(splits[2].c_str()));
					fiber.vertices.push_back(vrtx);
				}
			}
		}

		//printf("Yarn is initialized from the file. \n");
	}

	void Yarn::yarnCenter(const char *yarnfile, const char *yarnCenterfile) {

		std::ifstream fin;
		if (yarnfile != NULL)
			fin.open(yarnfile);

		//generate a one-ply yarn
		const int num_of_cores = omp_get_num_procs();
		Yarn yarn;
		const int ply_num = 1;
		int vrtx_num;
		std::string line;
		std::getline(fin, line);
		int fiber_num = atoi(line.c_str());
		yarn.plys.resize(ply_num);


		for (int p = 0; p < ply_num; ++p) {
			yarn.plys[p].fibers.resize(fiber_num);
			for (int f = 0; f < fiber_num; ++f) {
				Fiber &fiber = yarn.plys[p].fibers[f];
				fiber.clear(); //clear the vertices list 
				std::getline(fin, line);
				vrtx_num = atoi(line.c_str());
				for (int v = 0; v < vrtx_num; ++v) {
					std::getline(fin, line);
					std::vector<std::string> splits = split(line, ' ');
					vec3f vrtx(atof(splits[0].c_str()), atof(splits[1].c_str()), atof(splits[2].c_str()));
					fiber.vertices.push_back(vrtx);
				}
			}
		}
		//get the average
		FILE *fout;
		if (fopen_s(&fout, yarnCenterfile, "wt") == 0) {
			fprintf_s(fout, "%d \n", vrtx_num);
			for (int v = 0; v < vrtx_num; ++v) {
				vec3f sum_fibers = vec3f(0.f);
				for (int f = 0; f < fiber_num; ++f) {
					sum_fibers += yarn.plys[0].fibers[f].vertices[v];
				}
				sum_fibers /= static_cast<float>(fiber_num);
				fprintf_s(fout, "%.6f %.6f %.6f \n", sum_fibers.x, sum_fibers.y, sum_fibers.z);
				//fprintf_s(fout, "%.6f %.6f %.6f \n", 0.0, 0.0, sum_fibers.z);
			}
			fclose(fout);
		}

		printf("Centerfiber is written to the file. \n");
	}

	/* Treat each ply with the same settings */
	void Yarn::parse(const char *filename) {
		std::ifstream fin;
		if (filename != NULL)
			fin.open(filename);
		else
			fin.open(this->config_file.c_str());

		std::string line;
		while (std::getline(fin, line)) {
#ifdef VERBOSE
			std::cout << line << std::endl;
#endif
			std::vector<std::string> splits = split(line, ' ');
			if (splits.size() < 2)    continue;
			std::string p_name = splits[0];

			if (p_name == "ply_num:") {
				this->plys.resize(atoi(splits[1].c_str()));
			}
			else if (p_name == "fiber_num:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					int fiber_num = atoi(splits[1].c_str());
					this->plys[i].fibers.resize(fiber_num);
				}
#ifdef IMPROVED_FLYAWAYS
			}
			else if (p_name == "use_flyaways:") {
				this->use_flyaways = atoi(splits[1].c_str());
			}
			else if (p_name == "flyaway_hair_density:") {
				assert(this->plys.size());
				assert(splits.size() == 2);
				for (int i = 0; i < this->plys.size(); i++)
					this->plys[i].flyaway_hair_density = atof(splits[1].c_str());
			}
			else if (p_name == "flyaway_hair_ze:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < this->plys.size(); i++)
				{
					this->plys[i].flyaway_hair_ze_mu = atof(splits[1].c_str());
					this->plys[i].flyaway_hair_ze_sigma = atof(splits[2].c_str());
				}
			}
			else if (p_name == "flyaway_hair_r0:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < this->plys.size(); i++)
				{
					this->plys[i].flyaway_hair_r0_mu = atof(splits[1].c_str());
					this->plys[i].flyaway_hair_r0_sigma = atof(splits[2].c_str());
				}
			}
			else if (p_name == "flyaway_hair_re:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < this->plys.size(); i++)
				{
					this->plys[i].flyaway_hair_re_mu = atof(splits[1].c_str());
					this->plys[i].flyaway_hair_re_sigma = atof(splits[2].c_str());
				}
			}
			else if (p_name == "flyaway_hair_pe:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < this->plys.size(); i++)
				{
					this->plys[i].flyaway_hair_pe_mu = atof(splits[1].c_str());
					this->plys[i].flyaway_hair_pe_sigma = atof(splits[2].c_str());
				}
			}
			else if (p_name == "flyaway_loop_density:") {
				assert(this->plys.size());
				assert(splits.size() == 2);
				for (int i = 0; i < this->plys.size(); i++)
					this->plys[i].flyaway_loop_density = atof(splits[1].c_str());
			}
			else if (p_name == "flyaway_loop_r1:") {
				assert(this->plys.size());
				assert(splits.size() == 3);
				for (int i = 0; i < this->plys.size(); i++)
				{
					this->plys[i].flyaway_loop_r1_mu = atof(splits[1].c_str());
					this->plys[i].flyaway_loop_r1_sigma = atof(splits[2].c_str());
				}
#else
			}
			else if (p_name == "flyaway_num:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++)
					this->plys[i].flyaway_num = atoi(splits[1].c_str());
			}
			else if (p_name == "fly_step_size:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++)
					this->plys[i].fly_step_size = atof(splits[1].c_str());
#endif
			}
			else if (p_name == "z_step_size:") {
				this->z_step_size = atof(splits[1].c_str());
			}
			else if (p_name == "z_step_num:") {
				this->z_step_num = atof(splits[1].c_str());
			}
			else if (p_name == "yarn_clock_wise:") {
				this->clock_wise = atoi(splits[1].c_str());
			}
			else if (p_name == "fiber_clock_wise:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++)
					this->plys[i].clock_wise = atoi(splits[1].c_str());
			}
			else if (p_name == "use_migration:") {
				this->use_migration = atoi(splits[1].c_str());
			}
			else if (p_name == "yarn_alpha:") {
				this->yarn_alpha = atof(splits[1].c_str());
			}
			else if (p_name == "yarn_radius:") {
				this->yarn_radius = atof(splits[1].c_str());
			}
			else if (p_name == "epsilon:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].epsilon = atof(splits[1].c_str());
				}
			}
			else if (p_name == "R_max:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].R_max = atof(splits[1].c_str());
				}
			}
			else if (p_name == "beta:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].beta = atof(splits[1].c_str());
				}
			}
			else if (p_name == "alpha:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].alpha = atof(splits[1].c_str());
				}
			}
			else if (p_name == "s_i:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].s_i = atof(splits[1].c_str());
				}
			}
			else if (p_name == "rho_min:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].rho_min = atof(splits[1].c_str());
				}
			}
			else if (p_name == "rho_max:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].rho_max = atof(splits[1].c_str());
				}
			}
			else if (p_name == "ellipse_long:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].ellipse_long = atof(splits[1].c_str());
				}
			}
			else if (p_name == "ellipse_short:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].ellipse_short = atof(splits[1].c_str());
				}
#ifndef IMPROVED_FLYAWAYS
			}
			else if (p_name == "mu:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].mu = atof(splits[1].c_str());
				}
			}
			else if (p_name == "sigma:") {
				assert(this->plys.size());
				for (int i = 0; i < this->plys.size(); i++) {
					this->plys[i].sigma = atof(splits[1].c_str());
				}
#endif
			}
			else if (p_name == "z_curve_file:") {
				//this->z_curve.setFile(splits[1]); // Only one curve
			}
			else if (p_name == "z_curves_file:") {
				//this->z_curve.setFile(splits[1]);
				//this->z_curve.setMultiCurveFlag(true); // Multiple curves (e.g., woven)
			}
			else if (p_name == "aabb_min:") {
				std::string min_str = splits[1];
				std::vector<std::string> min_values = split(min_str.substr(1, min_str.size() - 2), ',');
				assert(min_values.size() == 3);
				for (int i = 0; i < 3; i++) {
					this->aabb_micro_ct.pMin[i] = atof(min_values[i].c_str());
				}
			}
			else if (p_name == "aabb_max:") {
				std::string max_str = splits[1];
				std::vector<std::string> max_values = split(max_str.substr(1, max_str.size() - 2), ',');
				assert(max_values.size() == 3);
				for (int i = 0; i < 3; i++) {
					this->aabb_micro_ct.pMax[i] = atof(max_values[i].c_str());
				}
			}
		}
		fin.close();
	}


	Yarn::Yarn() {}
	Yarn::~Yarn() {}

	void Yarn::simulate_ply() {


#define INDIVIDUAL_PLY 
		omp_init_lock(&this->lock);
		// Step1: Obtain center of yarn starting point 
#ifdef VERBOSE
		printf("Obtain center of yarn starting point...\n");
#endif
		const vec3f base_center = vec3f(0, 0, 0);
		const float base_radius = this->yarn_radius;
		this->aabb_procedural.reset();
		// Step2: Sample initial locations of ply-centers in normal plane around starting point
#ifdef VERBOSE
		printf("Sample initial locations of ply-centers in normal plane around starting point...\n");
#endif
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num;
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3f(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
		}

		// Step3: Sample initial fiber locations in normal plane around ply-centers using rejection sampling according to the distribution in Sec 4.1
#ifdef VERBOSE
		printf("Sample initial fiber locations in normal plane around ply-centers using rejection sampling according to the distribution in Sec 4.1...\n");
#endif

		const int num_of_cores = omp_get_num_procs();

		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
#pragma omp parallel for num_threads(num_of_cores) 
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				float radius = this->plys[i].sampleR();
				float theta = 2 * pi * (float)rand() / (RAND_MAX);
				float migration_theta = 2 * pi * (float)rand() / (RAND_MAX);
				fiber.init_radius = radius;
				fiber.init_theta = theta;
				fiber.init_migration_theta = migration_theta;
				fiber.init_vertex = this->plys[i].base_center +
					vec3f(radius * std::cosf(theta), radius * std::sinf(theta), 0);
			}

		}


		// Step4: Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center and rotating fiber positions around ply centers
#ifdef VERBOSE	
		printf("Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center and rotating fiber positions around ply centers...\n");
#endif


#ifdef IMPROVED_FLYAWAYS
		std::vector<std::vector<std::vector<float> > > rVals(ply_num);
#endif

#pragma omp parallel for num_threads(num_of_cores) 
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();

#ifdef IMPROVED_FLYAWAYS
			rVals[i].resize(fiber_num);
#else
			this->plys[i].flyaway_num = 0;
			this->plys[i].fly_fiber_num = 0;
#endif

			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				fiber.clear();

#ifdef PERTURB_FIBERS
				std::vector<float> perturbRatios;
				{
					const int nsteps = static_cast<int>(std::ceil(this->z_step_num));
					std::vector<int> eventLoc;
					for (int step_id = 0; step_id < nsteps; ++step_id)
						if ((float)rand() / (RAND_MAX) < PERTURB_FIBER_PROB)
							eventLoc.push_back(step_id);
					perturbRatios.resize(nsteps, 1.0f);
					if (!eventLoc.empty())
					{
						std::vector<int>::iterator it = eventLoc.begin();
						perturbRatios[*it] = 1.0f + PERTURB_FIBER_RATIO*((float)rand() / (RAND_MAX)- 0.5f);
						for (int j = 0; j < *it; ++j) perturbRatios[j] = perturbRatios[*it];
						while ((++it) != eventLoc.end())
						{
							perturbRatios[*it] = 1.0f + PERTURB_FIBER_RATIO*((float)rand() / (RAND_MAX)- 0.5f);
							float extent = static_cast<float>(*it - *(it - 1));
							for (int j = *(it - 1) + 1; j < *it; ++j)
							{
#if 0
								perturbRatios[j] = (perturbRatios[*(it - 1)] * (*it - j) + perturbRatios[*it] * (j - *(it - 1))) / extent;
#else
								float v = static_cast<float>(*it - j) / extent;
								v = std::sin(0.5f*pi*v);
								perturbRatios[j] = perturbRatios[*(it - 1)] * v + perturbRatios[*it] * (1.0f - v);
#endif
							}
						}
						for (int j = eventLoc.back() + 1; j < nsteps; ++j)
							perturbRatios[j] = perturbRatios[eventLoc.back()];
					}

					for (int j = 0; j < PERTURB_FIBER_SMOOTHING; ++j)
					{
						std::vector<float> perturbRatios0 = perturbRatios;
						for (int k = 1; k + 1 < nsteps; ++k)
							perturbRatios[k] = 0.25f*perturbRatios0[k - 1] + 0.5f*perturbRatios0[k] + 0.25f*perturbRatios0[k + 1];
					}
				}
#endif

#ifdef IMPROVED_FLYAWAYS
				rVals[i][f].clear();
#endif
				for (int step_id = 0; step_id < this->z_step_num; step_id++) {
					const float z = this->z_step_size * (step_id - this->z_step_num / 2.f);
					const float fiber_theta = this->plys[i].clock_wise ? -z * 2 * pi / this->plys[i].alpha : z * 2 * pi / this->plys[i].alpha;
					const float yarn_theta = this->clock_wise ? -z * 2 * pi / this->yarn_alpha : z * 2 * pi / this->yarn_alpha;
					float local_x, local_y, world_x, world_y;

					// Step5: Vary the distance of cross-sectional fiber positions to their ply center according to fiber migration Sec 4.2
					this->plys[i].helixXYZ(fiber.init_radius, fiber.init_theta, fiber_theta, use_migration, fiber.init_migration_theta, local_x, local_y);
#ifndef INDIVIDUAL_PLY
					// Step 6: Transform cross-sectional fiber positions according to strand compression Sec 4.3
					vec3f short_axis = nv::normalize(this->plys[i].base_center), long_axis = vec3f(-short_axis.y, short_axis.x, 0);
					vec3f local_p = vec3f(local_x, local_y, 0.f);
					float _local_x = nv::dot(local_p, short_axis), _local_y = nv::dot(local_p, long_axis);
					_local_x *= this->plys[i].ellipse_short;
					_local_y *= this->plys[i].ellipse_long;
					local_p = _local_x * short_axis + _local_y * long_axis;
					local_x = local_p.x;
					local_y = local_p.y;

#ifdef PERTURB_FIBERS
					local_x *= perturbRatios[step_id];
					local_y *= perturbRatios[step_id];
#endif

					float world_x_before_ply_rotation = local_x + this->plys[i].base_center.x;
					float world_y_before_ply_rotation = local_y + this->plys[i].base_center.y;
					world_x = world_x_before_ply_rotation * std::cosf(yarn_theta) - world_y_before_ply_rotation * std::sinf(yarn_theta);
					world_y = world_y_before_ply_rotation * std::cosf(yarn_theta) + world_x_before_ply_rotation * std::sinf(yarn_theta);
#else 
					const float balance_radius = std::sqrtf(this->plys[i].ellipse_short * this->plys[i].ellipse_long);
					local_x *= balance_radius;
					local_y *= balance_radius;
#ifdef IMPROVED_FLYAWAYS
					rVals[i][f].push_back(std::sqrt(local_x*local_x + local_y*local_y));
#endif
#ifdef PERTURB_FIBERS
					local_x *= perturbRatios[step_id];
					local_y *= perturbRatios[step_id];
#endif

					world_x = local_x;
					world_y = local_y;
#endif
					vec3f verIn = vec3f(world_x, world_y, z);
					fiber.vertices.push_back(verIn);
				}
			}
		}


#if 1 // 0 DISABLE FLYAWAY 

#ifdef IMPROVED_FLYAWAYS
		if (this->use_flyaways)
		{
			std::uniform_real_distribution<float> distrb1;
			std::normal_distribution<float> distrb2;
			std::mt19937 engine;
			engine.seed(/*1234*/rand());
#ifdef VERBOSE
			printf("Generating fly-away fibers...\n");
#endif
			const float sig_scale_hair = 0.75f, sig_scale_loop = 0.5f;
			const int min_loop_span = 10;

			float zextent = this->aabb_micro_ct.pMax.z - this->aabb_micro_ct.pMin.z;
			for (int i = 0; i < ply_num; ++i)
			{
				int nloop = static_cast<int>(std::floor(plys[i].flyaway_loop_density*zextent + 0.5f));
				if (nloop > 0)
				{

					std::vector<nv::vec2<int> > locs;
					int fiber_num = static_cast<int>(this->plys[i].fibers.size());
					for (int j = 0; j < fiber_num; ++j)
					{
						const Fiber &curFiber = this->plys[i].fibers[j];
						int totVtx = static_cast<int>(curFiber.vertices.size());
						for (int k = 1; k + 1 < totVtx; ++k)
							if (rVals[i][j][k] > rVals[i][j][k - 1] && rVals[i][j][k] > rVals[i][j][k + 1])
								locs.push_back(nv::vec2<int>(j, k));
					}
					std::random_shuffle(locs.begin(), locs.end());

					for (int j = 0; j < nloop && j < static_cast<int>(locs.size()); ++j)
					{
						int fid = locs[j][0];
						const std::vector<float> curRs = rVals[i][fid];
						Fiber &curFiber = this->plys[i].fibers[fid];
						int totVtx = static_cast<int>(curFiber.vertices.size());

						const int k = locs[j][1];
						int k0 = std::max(k - min_loop_span, 0);
						int k1 = std::min(k + min_loop_span, totVtx - 1);
						while (k0 > 0 && curRs[k0 - 1] < curRs[k0]) --k0;
						while (k1 + 1 < totVtx && curRs[k1 + 1] < curRs[k1]) ++k1;

						float r1;
						for (; ; )
						{
							r1 = plys[i].flyaway_loop_r1_mu + sig_scale_loop*plys[i].flyaway_loop_r1_sigma*distrb2(engine);
							if (r1 > 1.05f*curRs[k]) break;
						}

						float ratio = r1 / curRs[k];
						for (int t = k0 + 1; t <= k; ++t)
						{
							//float v = 1.0f + (ratio - 1.0f)*static_cast<float>(t - k0)/static_cast<float>(k - k0);
							float v = 1.0f + (ratio - 1.0f)*std::sin(0.5f*pi*static_cast<float>(t - k0) / static_cast<float>(k - k0));
							curFiber.vertices[t].x *= v;
							curFiber.vertices[t].y *= v;
						}
						for (int t = k1 - 1; t > k; --t)
						{
							//float v = 1.0f + (ratio - 1.0f)*static_cast<float>(t - k1)/static_cast<float>(k - k1);
							float v = 1.0f + (ratio - 1.0f)*std::sin(0.5f*pi*static_cast<float>(t - k1) / static_cast<float>(k - k1));
							curFiber.vertices[t].x *= v;
							curFiber.vertices[t].y *= v;
						}
					}
				}

				int nhair = static_cast<int>(std::floor(plys[i].flyaway_hair_density*zextent + 0.5f));
				for (int j = 0; j < nhair; )
				{
					Fiber fiber;
					float z0 = this->aabb_micro_ct.pMin.z + distrb1(engine)*zextent;
					float ze = plys[i].flyaway_hair_ze_mu + sig_scale_hair*plys[i].flyaway_hair_ze_sigma*distrb2(engine);
					float r0 = plys[i].flyaway_hair_r0_mu + sig_scale_hair*plys[i].flyaway_hair_r0_sigma*distrb2(engine);
					float re = plys[i].flyaway_hair_re_mu + sig_scale_hair*plys[i].flyaway_hair_re_sigma*distrb2(engine);
					float p0 = 2.0f*pi*distrb1(engine);
					float pe = plys[i].flyaway_hair_pe_mu + sig_scale_hair*plys[i].flyaway_hair_pe_sigma*distrb2(engine);

					/* Extrapolation */

					float r0_e = 0.0f, re_e = r0 + re;
					float z0_e = z0 - ze*r0 / re, ze_e = ze + ze*r0 / re;
					float p0_e = p0 - pe*r0 / re, pe_e = pe + pe*r0 / re;

					int nstep = 100;
					std::vector<vec3f> vars;
					for (int k = 0; k <= nstep; ++k)
					{
						vec3f cur;
						cur[0] = r0_e + re_e*static_cast<float>(k) / static_cast<float>(nstep);
						cur[1] = z0_e + ze_e*static_cast<float>(k) / static_cast<float>(nstep);
						cur[2] = p0_e + pe_e*static_cast<float>(k) / static_cast<float>(nstep);
						vars.push_back(cur);
					}

#ifdef PERTURB_FIBERS
					/* Perturb parameters */
					std::vector<float> perturbRatios;
					for (int k = 0; k < 1; ++k)
					{
						std::vector<int> eventLoc;
						for (int step_id = 0; step_id <= nstep; ++step_id)
							if ((float)rand() / (RAND_MAX)< 0.2f /* PERTURB_FIBER_PROB */)
								eventLoc.push_back(step_id);
						perturbRatios.resize(nstep + 1, 1.0f);
						if (!eventLoc.empty())
						{
							std::vector<int>::iterator it = eventLoc.begin();
							perturbRatios[*it] = 1.0f + 0.1f /* PERTURB_FIBER_RATIO */ *((float)rand() / (RAND_MAX) - 0.5f);
							for (int t = 0; t < *it; ++t) perturbRatios[t] = perturbRatios[*it];
							while ((++it) != eventLoc.end())
							{
								perturbRatios[*it] = 1.0f + 0.1f /* PERTURB_FIBER_RATIO */ *((float)rand() / (RAND_MAX) - 0.5f);
								float extent = static_cast<float>(*it - *(it - 1));
								for (int t = *(it - 1) + 1; t < *it; ++t)
								{
#if 0
									perturbRatios[t] = (perturbRatios[*(it - 1)] * (*it - t) + perturbRatios[*it] * (t - *(it - 1))) / extent;
#else
									float v = static_cast<float>(*it - t) / extent;
									v = std::sin(0.5f*pi*v);
									perturbRatios[t] = perturbRatios[*(it - 1)] * v + perturbRatios[*it] * (1.0f - v);
#endif
								}
							}
							for (int t = eventLoc.back() + 1; t <= nstep; ++t)
								perturbRatios[t] = perturbRatios[eventLoc.back()];
						}

						for (int t = 0; t < PERTURB_FIBER_SMOOTHING; ++t)
						{
							std::vector<float> perturbRatios0 = perturbRatios;
							for (int o = 1; o + 1 <= nstep; ++o)
								perturbRatios[o] = 0.25f*perturbRatios0[o - 1] + 0.5f*perturbRatios0[o] + 0.25f*perturbRatios0[o + 1];
						}

						for (int t = 0; t <= nstep; ++t) vars[t][k] *= perturbRatios[t];
					}
#endif	
					/* Creating fiber curve */

					for (int k = 0; k <= nstep; ++k)
					{
						const vec3f &cur = vars[k];
						vec3f pos;
						pos[0] = cur[0] * std::cos(cur[2]);
						pos[1] = cur[0] * std::sin(cur[2]);
						pos[2] = cur[1];
						// Crop flyaway fibers using the ply's bounding box
						if (pos[2] < this->aabb_micro_ct.pMin.z || pos[2] > this->aabb_micro_ct.pMax.z)
							break;
						fiber.vertices.push_back(pos);
					}
					if (fiber.vertices.size() > 1)
					{
						this->plys[i].fibers.push_back(fiber);
						++j;
					}
				}
#ifdef VERBOSE
				printf("    Ply #%d: %d type-hair fibers, %d type-loop fibers\n", i, nhair, nloop);
#endif
			}
		}
#endif

#endif // DISABLE FLYAWAY 



		for (int i = 0; i < ply_num; ++i)
			for (auto it = this->plys[i].fibers.begin(); it != this->plys[i].fibers.end(); ++it)
				for (auto it2 = it->vertices.begin(); it2 != it->vertices.end(); ++it2) {
					this->aabb_procedural.grow(*it2);
				}


		omp_destroy_lock(&this->lock);




		std::cout << "Checking..." << std::endl;
		int bad_count = 0;
		for (int i = 0; i < this->plys.size(); i++)
		{
			const int fiberNum = this->plys[i].fibers.size();
			for (int f = 0; f < fiberNum; f++)
			{
				const int vertexNum = this->plys[i].fibers[f].vertices.size();
				for (int v = 1; v < vertexNum - 1; v++)
				{
					vec3f prev = this->plys[i].fibers[f].vertices[v - 1];
					vec3f curr = this->plys[i].fibers[f].vertices[v];
					vec3f next = this->plys[i].fibers[f].vertices[v + 1];
					vec3f dir1 = nv::normalize(next - curr);
					vec3f dir2 = nv::normalize(curr - prev);
					if (nv::dot(dir1, dir2) < 0.5) {
						bad_count++;
					}
				}
			}
		}
		std::cout << "Bad count = " << bad_count << std::endl;

	}

#if 0
	void Yarn::simulate_ply() {

#define INDIVIDUAL_PLY 
		omp_init_lock(&this->lock);
		// Step1: Obtain center of yarn starting point 
#ifdef VERBOSE
		printf("Obtain center of yarn starting point...\n");
#endif
		const vec3f base_center = vec3f(0, 0, 0);
		const float base_radius = this->yarn_radius;
		this->aabb_procedural.reset();
		// Step2: Sample initial locations of ply-centers in normal plane around starting point
#ifdef VERBOSE
		printf("Sample initial locations of ply-centers in normal plane around starting point...\n");
#endif
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num;
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3f(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
		}

		// Step3: Sample initial fiber locations in normal plane around ply-centers using rejection sampling according to the distribution in Sec 4.1
#ifdef VERBOSE
		printf("Sample initial fiber locations in normal plane around ply-centers using rejection sampling according to the distribution in Sec 4.1...\n");
#endif

		const int num_of_cores = omp_get_num_procs();

		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
//#pragma omp parallel for num_threads(num_of_cores) 
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];

				float radius;
//#pragma omp critical
				radius = this->plys[i].sampleR();

				float theta = 2 * pi * (float)rand() / (RAND_MAX);
				float migration_theta = 2 * pi * (float)rand() / (RAND_MAX);
				fiber.init_radius = radius;
				fiber.init_theta = theta;
				fiber.init_migration_theta = migration_theta;
				fiber.init_vertex = this->plys[i].base_center +
					vec3f(radius * std::cosf(theta), radius * std::sinf(theta), 0);
			}

		}


		// Step4: Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center and rotating fiber positions around ply centers
#ifdef VERBOSE	
		printf("Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center and rotating fiber positions around ply centers...\n");
#endif


#ifdef IMPROVED_FLYAWAYS
		std::vector<std::vector<std::vector<float> > > rVals(ply_num);
#endif

//#pragma omp parallel for num_threads(num_of_cores) 
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();

#ifdef IMPROVED_FLYAWAYS
			rVals[i].resize(fiber_num);
#else
			this->plys[i].flyaway_num = 0;
			this->plys[i].fly_fiber_num = 0;
#endif

			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				fiber.clear();

#ifdef PERTURB_FIBERS
				std::vector<float> perturbRatios;
				{
					const int nsteps = static_cast<int>(std::ceil(this->z_step_num));
					std::vector<int> eventLoc;
					for (int step_id = 0; step_id < nsteps; ++step_id)
						if ((float)rand() / (RAND_MAX) < PERTURB_FIBER_PROB)
							eventLoc.push_back(step_id);
					perturbRatios.resize(nsteps, 1.0f);
					if (!eventLoc.empty())
					{
						std::vector<int>::iterator it = eventLoc.begin();
						perturbRatios[*it] = 1.0f + PERTURB_FIBER_RATIO*((float)rand() / (RAND_MAX)- 0.5f);
						for (int j = 0; j < *it; ++j) perturbRatios[j] = perturbRatios[*it];
						while ((++it) != eventLoc.end())
						{
							perturbRatios[*it] = 1.0f + PERTURB_FIBER_RATIO*((float)rand() / (RAND_MAX)- 0.5f);
							float extent = static_cast<float>(*it - *(it - 1));
							for (int j = *(it - 1) + 1; j < *it; ++j)
							{
#if 0
								perturbRatios[j] = (perturbRatios[*(it - 1)] * (*it - j) + perturbRatios[*it] * (j - *(it - 1))) / extent;
#else
								float v = static_cast<float>(*it - j) / extent;
								v = std::sin(0.5f*pi*v);
								perturbRatios[j] = perturbRatios[*(it - 1)] * v + perturbRatios[*it] * (1.0f - v);
#endif
							}
						}
						for (int j = eventLoc.back() + 1; j < nsteps; ++j)
							perturbRatios[j] = perturbRatios[eventLoc.back()];
					}

					for (int j = 0; j < PERTURB_FIBER_SMOOTHING; ++j)
					{
						std::vector<float> perturbRatios0 = perturbRatios;
						for (int k = 1; k + 1 < nsteps; ++k)
							perturbRatios[k] = 0.25f*perturbRatios0[k - 1] + 0.5f*perturbRatios0[k] + 0.25f*perturbRatios0[k + 1];
					}
				}
#endif

#ifdef IMPROVED_FLYAWAYS
				rVals[i][f].clear();
#endif
				for (int step_id = 0; step_id < this->z_step_num; step_id++) {
					const float z = this->z_step_size * (step_id - this->z_step_num / 2.f);
					const float fiber_theta = this->plys[i].clock_wise ? -z * 2 * pi / this->plys[i].alpha : z * 2 * pi / this->plys[i].alpha;
					const float yarn_theta = this->clock_wise ? -z * 2 * pi / this->yarn_alpha : z * 2 * pi / this->yarn_alpha;
					float local_x, local_y, world_x, world_y;

					// Step5: Vary the distance of cross-sectional fiber positions to their ply center according to fiber migration Sec 4.2
					this->plys[i].helixXYZ(fiber.init_radius, fiber.init_theta, fiber_theta, use_migration, fiber.init_migration_theta, local_x, local_y);
#ifndef INDIVIDUAL_PLY
					// Step 6: Transform cross-sectional fiber positions according to strand compression Sec 4.3
					vec3f short_axis = nv::normalize(this->plys[i].base_center), long_axis = vec3f(-short_axis.y, short_axis.x, 0);
					vec3f local_p = vec3f(local_x, local_y, 0.f);
					float _local_x = nv::dot(local_p, short_axis), _local_y = nv::dot(local_p, long_axis);
					_local_x *= this->plys[i].ellipse_short;
					_local_y *= this->plys[i].ellipse_long;
					local_p = _local_x * short_axis + _local_y * long_axis;
					local_x = local_p.x;
					local_y = local_p.y;

#ifdef PERTURB_FIBERS
					local_x *= perturbRatios[step_id];
					local_y *= perturbRatios[step_id];
#endif

					float world_x_before_ply_rotation = local_x + this->plys[i].base_center.x;
					float world_y_before_ply_rotation = local_y + this->plys[i].base_center.y;
					world_x = world_x_before_ply_rotation * std::cosf(yarn_theta) - world_y_before_ply_rotation * std::sinf(yarn_theta);
					world_y = world_y_before_ply_rotation * std::cosf(yarn_theta) + world_x_before_ply_rotation * std::sinf(yarn_theta);
#else 
					const float balance_radius = std::sqrtf(this->plys[i].ellipse_short * this->plys[i].ellipse_long);
					local_x *= balance_radius;
					local_y *= balance_radius;
#ifdef IMPROVED_FLYAWAYS
					rVals[i][f].push_back(std::sqrt(local_x*local_x + local_y*local_y));
#endif
#ifdef PERTURB_FIBERS
					local_x *= perturbRatios[step_id];
					local_y *= perturbRatios[step_id];
#endif

					world_x = local_x;
					world_y = local_y;
#endif
					vec3f verIn = vec3f(world_x, world_y, z);
					fiber.vertices.push_back(verIn);
				}
			}
		}


#if 1 // 0 DISABLE FLYAWAY 

#ifdef IMPROVED_FLYAWAYS
		if (this->use_flyaways)
		{
			std::uniform_real_distribution<float> distrb1;
			std::normal_distribution<float> distrb2;
			std::mt19937 engine;
			engine.seed(/*1234*/rand());
#ifdef VERBOSE
			printf("Generating fly-away fibers...\n");
#endif
			const float sig_scale_hair = 0.75f, sig_scale_loop = 0.5f;
			const int min_loop_span = 10;

			float zextent = this->aabb_micro_ct.pMax.z - this->aabb_micro_ct.pMin.z;
			for (int i = 0; i < ply_num; ++i)
			{
				int nloop = static_cast<int>(std::floor(plys[i].flyaway_loop_density*zextent + 0.5f));
				if (nloop > 0)
				{

					std::vector<nv::vec2<int> > locs;
					int fiber_num = static_cast<int>(this->plys[i].fibers.size());
					for (int j = 0; j < fiber_num; ++j)
					{
						const Fiber &curFiber = this->plys[i].fibers[j];
						int totVtx = static_cast<int>(curFiber.vertices.size());
						for (int k = 1; k + 1 < totVtx; ++k)
							if (rVals[i][j][k] > rVals[i][j][k - 1] && rVals[i][j][k] > rVals[i][j][k + 1])
								locs.push_back(nv::vec2<int>(j, k));
					}
					std::random_shuffle(locs.begin(), locs.end());

					for (int j = 0; j < nloop && j < static_cast<int>(locs.size()); ++j)
					{
						int fid = locs[j][0];
						const std::vector<float> curRs = rVals[i][fid];
						Fiber &curFiber = this->plys[i].fibers[fid];
						int totVtx = static_cast<int>(curFiber.vertices.size());

						const int k = locs[j][1];
						int k0 = std::max(k - min_loop_span, 0);
						int k1 = std::min(k + min_loop_span, totVtx - 1);
						while (k0 > 0 && curRs[k0 - 1] < curRs[k0]) --k0;
						while (k1 + 1 < totVtx && curRs[k1 + 1] < curRs[k1]) ++k1;

						float r1;
						for (; ; )
						{
							r1 = plys[i].flyaway_loop_r1_mu + sig_scale_loop*plys[i].flyaway_loop_r1_sigma*distrb2(engine);
							if (r1 > 1.05f*curRs[k]) break;
						}

						float ratio = r1 / curRs[k];
						for (int t = k0 + 1; t <= k; ++t)
						{
							//float v = 1.0f + (ratio - 1.0f)*static_cast<float>(t - k0)/static_cast<float>(k - k0);
							float v = 1.0f + (ratio - 1.0f)*std::sin(0.5f*pi*static_cast<float>(t - k0) / static_cast<float>(k - k0));
							curFiber.vertices[t].x *= v;
							curFiber.vertices[t].y *= v;
						}
						for (int t = k1 - 1; t > k; --t)
						{
							//float v = 1.0f + (ratio - 1.0f)*static_cast<float>(t - k1)/static_cast<float>(k - k1);
							float v = 1.0f + (ratio - 1.0f)*std::sin(0.5f*pi*static_cast<float>(t - k1) / static_cast<float>(k - k1));
							curFiber.vertices[t].x *= v;
							curFiber.vertices[t].y *= v;
						}
					}
				}

				int nhair = static_cast<int>(std::floor(plys[i].flyaway_hair_density*zextent + 0.5f));
				for (int j = 0; j < nhair; )
				{
					Fiber fiber;
					float z0 = this->aabb_micro_ct.pMin.z + distrb1(engine)*zextent;
					float ze = plys[i].flyaway_hair_ze_mu + sig_scale_hair*plys[i].flyaway_hair_ze_sigma*distrb2(engine);
					float r0 = plys[i].flyaway_hair_r0_mu + sig_scale_hair*plys[i].flyaway_hair_r0_sigma*distrb2(engine);
					float re = plys[i].flyaway_hair_re_mu + sig_scale_hair*plys[i].flyaway_hair_re_sigma*distrb2(engine);
					float p0 = 2.0f*pi*distrb1(engine);
					float pe = plys[i].flyaway_hair_pe_mu + sig_scale_hair*plys[i].flyaway_hair_pe_sigma*distrb2(engine);

					/* Extrapolation */

					float r0_e = 0.0f, re_e = r0 + re;
					float z0_e = z0 - ze*r0 / re, ze_e = ze + ze*r0 / re;
					float p0_e = p0 - pe*r0 / re, pe_e = pe + pe*r0 / re;

					int nstep = 100;
					std::vector<vec3f> vars;
					for (int k = 0; k <= nstep; ++k)
					{
						vec3f cur;
						cur[0] = r0_e + re_e*static_cast<float>(k) / static_cast<float>(nstep);
						cur[1] = z0_e + ze_e*static_cast<float>(k) / static_cast<float>(nstep);
						cur[2] = p0_e + pe_e*static_cast<float>(k) / static_cast<float>(nstep);
						vars.push_back(cur);
					}

#ifdef PERTURB_FIBERS
					/* Perturb parameters */
					std::vector<float> perturbRatios;
					for (int k = 0; k < 1; ++k)
					{
						std::vector<int> eventLoc;
						for (int step_id = 0; step_id <= nstep; ++step_id)
							if ((float)rand() / (RAND_MAX) < 0.2f /* PERTURB_FIBER_PROB */)
								eventLoc.push_back(step_id);
						perturbRatios.resize(nstep + 1, 1.0f);
						if (!eventLoc.empty())
						{
							std::vector<int>::iterator it = eventLoc.begin();
							perturbRatios[*it] = 1.0f + 0.1f /* PERTURB_FIBER_RATIO */ *((float)rand() / (RAND_MAX)- 0.5f);
							for (int t = 0; t < *it; ++t) perturbRatios[t] = perturbRatios[*it];
							while ((++it) != eventLoc.end())
							{
								perturbRatios[*it] = 1.0f + 0.1f /* PERTURB_FIBER_RATIO */ *((float)rand() / (RAND_MAX)- 0.5f);
								float extent = static_cast<float>(*it - *(it - 1));
								for (int t = *(it - 1) + 1; t < *it; ++t)
								{
#if 0
									perturbRatios[t] = (perturbRatios[*(it - 1)] * (*it - t) + perturbRatios[*it] * (t - *(it - 1))) / extent;
#else
									float v = static_cast<float>(*it - t) / extent;
									v = std::sin(0.5f*pi*v);
									perturbRatios[t] = perturbRatios[*(it - 1)] * v + perturbRatios[*it] * (1.0f - v);
#endif
								}
							}
							for (int t = eventLoc.back() + 1; t <= nstep; ++t)
								perturbRatios[t] = perturbRatios[eventLoc.back()];
						}

						for (int t = 0; t < PERTURB_FIBER_SMOOTHING; ++t)
						{
							std::vector<float> perturbRatios0 = perturbRatios;
							for (int o = 1; o + 1 <= nstep; ++o)
								perturbRatios[o] = 0.25f*perturbRatios0[o - 1] + 0.5f*perturbRatios0[o] + 0.25f*perturbRatios0[o + 1];
						}

						for (int t = 0; t <= nstep; ++t) vars[t][k] *= perturbRatios[t];
					}
#endif	
					/* Creating fiber curve */

					for (int k = 0; k <= nstep; ++k)
					{
						const vec3f &cur = vars[k];
						vec3f pos;
						pos[0] = cur[0] * std::cos(cur[2]);
						pos[1] = cur[0] * std::sin(cur[2]);
						pos[2] = cur[1];
						// Crop flyaway fibers using the ply's bounding box
						if (pos[2] < this->aabb_micro_ct.pMin.z || pos[2] > this->aabb_micro_ct.pMax.z)
							break;
						fiber.vertices.push_back(pos);
					}
					if (fiber.vertices.size() > 1)
					{
						this->plys[i].fibers.push_back(fiber);
						++j;
					}
				}
#ifdef VERBOSE
				printf("    Ply #%d: %d type-hair fibers, %d type-loop fibers\n", i, nhair, nloop);
#endif
			}
		}
#endif

#endif // DISABLE FLYAWAY 



		for (int i = 0; i < ply_num; ++i)
			for (auto it = this->plys[i].fibers.begin(); it != this->plys[i].fibers.end(); ++it)
				for (auto it2 = it->vertices.begin(); it2 != it->vertices.end(); ++it2) {
					this->aabb_procedural.grow(*it2);
				}


		omp_destroy_lock(&this->lock);




		std::cout << "Checking..." << std::endl;
		int bad_count = 0;
		for (int i = 0; i < this->plys.size(); i++)
		{
			const int fiberNum = this->plys[i].fibers.size();
			for (int f = 0; f < fiberNum; f++)
			{
				const int vertexNum = this->plys[i].fibers[f].vertices.size();
				for (int v = 1; v < vertexNum - 1; v++)
				{
					vec3f prev = this->plys[i].fibers[f].vertices[v - 1];
					vec3f curr = this->plys[i].fibers[f].vertices[v];
					vec3f next = this->plys[i].fibers[f].vertices[v + 1];
					vec3f dir1 = nv::normalize(next - curr);
					vec3f dir2 = nv::normalize(curr - prev);
					if (nv::dot(dir1, dir2) < 0.5) {
						bad_count++;
					}
				}
			}
		}
		std::cout << "Bad count = " << bad_count << std::endl;

	}
#endif
	void Yarn::roll_plys(const int K, const std::string &ply_fn, const std::string &fiber_fn) {
		const int num_of_cores = omp_get_num_procs();
#ifdef VERBOSE
		std::cout << "Rolling plys into yarn... " << std::endl;
#endif
		typedef std::vector<vec3f> fiber_t;
		typedef std::vector<fiber_t> ply_t;

		std::vector<ply_t> yarn(K);
		int total_fiber_num = 0; std::vector<int> ply_fiber_num(K);

		for (int i = 0; i < K; i++) {
			std::string filename = ply_fn.substr(0, ply_fn.find('.')) + std::to_string((long long)i) + ".txt";
#ifdef VERBOSE
			std::cout << "Loading unrolled ply = " << filename << std::endl;
#endif
			std::ifstream fin(filename.c_str());
			int fiberNum; fin >> fiberNum;
			ply_fiber_num[i] = fiberNum;
			total_fiber_num += fiberNum;
			yarn[i].resize(fiberNum);
			for (int f = 0; f < fiberNum; f++) {
				int vertexNum; fin >> vertexNum;
				yarn[i][f].resize(vertexNum);
				for (int v = 0; v < vertexNum; v++) {
					fin >> yarn[i][f][v].x >> yarn[i][f][v].y >> yarn[i][f][v].z;
				}
			}
			fin.close();
		}

#ifdef VERBOSE
		printf("Obtain center of yarn starting point...\n");
#endif
		const vec3f base_center = vec3f(0, 0, 0);
		const float base_radius = this->yarn_radius;

#ifdef VERBOSE
		printf("Sample initial locations of ply-centers in normal plane around starting point...\n");
#endif
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num;
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3f(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
		}

		const float yarn_alpha = this->yarn_alpha;
		const int   yarn_clock = this->clock_wise;

#ifdef VERBOSE
		printf("Follow cross-section vertices along yarn center paths, while rotating ply centers around the yarn center...\n");
#endif
#pragma omp parallel for num_threads(num_of_cores) 
		for (int i = 0; i < ply_num; i++) {

			const float e_l = this->plys[i].ellipse_long, e_s = this->plys[i].ellipse_short, b_r = sqrtf(e_l * e_s);

			const int fiberNum = yarn[i].size();
			for (int f = 0; f < fiberNum; f++) {
				const int vertexNum = yarn[i][f].size();
				for (int v = 0; v < vertexNum; v++) {
					const float z = yarn[i][f][v].z;
					const float yarn_theta = yarn_clock ? -z * 2 * pi / yarn_alpha : z * 2 * pi / yarn_alpha;
					float local_x, local_y, world_x, world_y;

					local_x = yarn[i][f][v].x / b_r;
					local_y = yarn[i][f][v].y / b_r;


					vec3f short_axis = nv::normalize(this->plys[i].base_center), long_axis = vec3f(-short_axis.y, short_axis.x, 0);
					vec3f local_p = vec3f(local_x, local_y, 0.f);
					float _local_x = nv::dot(local_p, short_axis), _local_y = nv::dot(local_p, long_axis);
					_local_x *= this->plys[i].ellipse_short;
					_local_y *= this->plys[i].ellipse_long;
					local_p = _local_x * short_axis + _local_y * long_axis;
					local_x = local_p.x;
					local_y = local_p.y;
					float world_x_before_ply_rotation = local_x + this->plys[i].base_center.x;
					float world_y_before_ply_rotation = local_y + this->plys[i].base_center.y;
					world_x = world_x_before_ply_rotation * std::cosf(yarn_theta) - world_y_before_ply_rotation * std::sinf(yarn_theta);
					world_y = world_y_before_ply_rotation * std::cosf(yarn_theta) + world_x_before_ply_rotation * std::sinf(yarn_theta);

					yarn[i][f][v].x = world_x;
					yarn[i][f][v].y = world_y;
				}
			}

		}

#ifdef VERBOSE
		std::cout << "Writing final yarn file... " << std::endl;
#endif
		std::ofstream fout(fiber_fn.c_str());

		fout << total_fiber_num << std::endl;

		for (int i = 0; i < K; i++) {
			const int fiberNum = yarn[i].size();
			for (int f = 0; f < fiberNum; f++) {
				int vertexNum = yarn[i][f].size();
				fout << vertexNum << std::endl;
				for (int v = 0; v < vertexNum; v++) {
					fout << yarn[i][f][v].x << " " << yarn[i][f][v].y << " " << yarn[i][f][v].z << std::endl;
				}
			}
		}


		fout.close();
#ifdef VERBOSE
		std::cout << "Done!" << std::endl;
#endif
	}


	void Yarn::yarn_simulate() {

		std::cout << "step1: yarn center ...\n";
		const float base_radius = this->yarn_radius;
		this->aabb_procedural.reset();

		std::cout << "step2: initial plys centers ... \n";
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num; // place ply centers in a circle with radius yarn-radius/2 
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3f(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
		}

		const int num_of_cores = omp_get_num_procs();

		std::cout << "step3: initial plys locations ... \n";
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
//#pragma omp parallel for num_threads(num_of_cores) 
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				
				float radius;
//#pragma omp critical
				radius = this->plys[i].sampleR();

				float theta = 2 * pi * (float)rand() / (RAND_MAX);
				float migration_theta = 2 * pi * (float)rand() / (RAND_MAX);
				fiber.init_radius = radius;
				fiber.init_theta = theta;
				fiber.init_migration_theta = migration_theta;
				fiber.init_vertex = this->plys[i].base_center +
					vec3f(radius * std::cosf(theta), radius * std::sinf(theta), 0.2);
			}
		}

		std::cout << "step4-5-6: rotate ply-centers around yarn-center and fibers around ply-centers and apply the compression ... \n";
//#pragma omp parallel for num_threads(num_of_cores)
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();

			// generate all fibers around ply-center
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				fiber.clear(); //clear the vertices list 
				for (int step_id = 0; step_id < this->z_step_num; step_id++) {
					const float z = this->z_step_size * (step_id - this->z_step_num / 2.f); // devided by 2 Bcuz yarn lies between neg and pos z
					const float fiber_theta = this->plys[i].clock_wise ? -z * 2 * pi / this->plys[i].alpha : z * 2 * pi / this->plys[i].alpha;
					const float yarn_theta = this->clock_wise ? -z * 2 * pi / this->yarn_alpha : z * 2 * pi / this->yarn_alpha;
					float local_x, local_y, world_x, world_y;

					// std::cout << "step5: generate ply after fiber migration ... \n";
					// translate positions to rotate around ply center (use fiber_theta)
					this->plys[i].helixXYZ(fiber.init_radius, fiber.init_theta, fiber_theta, use_migration, fiber.init_migration_theta, local_x, local_y);

					// std::cout << "step5: eliptical compression ... \n";
					vec3f short_axis = nv::normalize(this->plys[i].base_center), long_axis = vec3f(-short_axis.y, short_axis.x, 0);
					// short axis is [ply-center - 0] because its a vector pointing to origin starting at first cross section, the long axis is perpendicular to this.
					vec3f local_p = vec3f(local_x, local_y, 0.f);
					float _local_x = nv::dot(local_p, short_axis), _local_y = nv::dot(local_p, long_axis);
					// scale the shape of cross section
					_local_x *= this->plys[i].ellipse_short;
					_local_y *= this->plys[i].ellipse_long;
					local_p = _local_x * short_axis + _local_y * long_axis;
					local_x = local_p.x;
					local_y = local_p.y;

					// 1. translate positions to rotate around yarn center (use yarn_theta)
					float world_x_before_ply_rotation = local_x + this->plys[i].base_center.x;
					float world_y_before_ply_rotation = local_y + this->plys[i].base_center.y;

					// 2. rotate positions [x' y'] = [cos sin; -sin cos][x y]
					world_x = world_x_before_ply_rotation * std::cosf(yarn_theta) - world_y_before_ply_rotation * std::sinf(yarn_theta);
					world_y = world_y_before_ply_rotation * std::cosf(yarn_theta) + world_x_before_ply_rotation * std::sinf(yarn_theta);

					vec3f verIn = vec3f(world_x, world_y, z), verOut;
					verOut = verIn;

					this->aabb_procedural.grow(verOut);
					if (this->aabb_micro_ct.in(verOut))
						fiber.vertices.push_back(verOut);
				}
			}
		}
	} // yarn_simulate

	void Yarn::readCompressFile_A(const char* compress_S, std::vector<Eigen::Matrix2f> &all_A) {

		std::ifstream fin;
		fin.open(compress_S);
		std::string line;
		//std::getline(fin, line);

		int i = 0;
		while (std::getline(fin, line)) {
			std::vector<std::string> splits = split(line, ' ');
			float S00 = atof(splits[0].c_str());
			float S01 = atof(splits[1].c_str());
			float S10 = atof(splits[2].c_str());
			float S11 = atof(splits[3].c_str());
			const int limit = 1000;
			Eigen::Matrix2f A;
			if (S00 > limit || S01 > limit || S10 > limit || S11 > limit)
				A = Eigen::Matrix2f::Identity();
			else
				A << S00, S01, S10, S11;

			all_A.push_back(A);
			i++;
		}
	}


	void Yarn::flip_yarn(const char* global_rot) {
		std::cout << "step7: rotate yarn cross-sections ..." << std::endl;

		double zMin = std::numeric_limits<double>::max(), zMax = std::numeric_limits<double>::lowest();
		for (const auto &ply : plys)
			for (const auto &fiber : ply.fibers)
				for (const auto &vertex : fiber.vertices) {
					zMin = std::min(zMin, static_cast<double>(vertex.z));
					zMax = std::max(zMax, static_cast<double>(vertex.z));
				}
		double zSpan = zMax - zMin;

		const int ply_num = this->plys.size();
		const int num_of_cores = omp_get_num_procs();

		std::vector<Eigen::Matrix2f> R;
		readCompressFile_A(global_rot, R);

		// change the yarn cross-sections
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
#pragma omp parallel for num_threads(num_of_cores) 
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				const int vertices_num = this->plys[i].fibers[f].vertices.size();

				for (int v = 0; v < vertices_num; v++) {

					fiber.vertices[v].x = -1.0 * fiber.vertices[v].x;
					fiber.vertices[v].y = 1.0 * fiber.vertices[v].y;
					fiber.vertices[v].z = 1.0 * fiber.vertices[v].z;
				}
			}
		}
	} // rotate_yarn


	void Yarn::compress_yarn_A(const char* compress_S, const char* global_rot) {
		std::cout << "step7: compress yarn cross-sections ..." << std::endl;

		double zMin = std::numeric_limits<double>::max(), zMax = std::numeric_limits<double>::lowest();
		for (const auto &ply : plys)
			for (const auto &fiber : ply.fibers)
				for (const auto &vertex : fiber.vertices) {
					zMin = std::min(zMin, static_cast<double>(vertex.z));
					zMax = std::max(zMax, static_cast<double>(vertex.z));
				}
		double zSpan = zMax - zMin;
		//double curveLength = centerCurve.totalLength();

		const int ply_num = this->plys.size();

		std::vector<Eigen::Matrix2f> A;
		std::vector<Eigen::Matrix2f> R;
		readCompressFile_A(compress_S, A); //stretching factor S
		if (global_rot != "") {
			readCompressFile_A(global_rot, R); //rotation factor R
		}



		if (A.size() != this->z_step_num)
			std::cout << "# compress params: " << A.size() << ", # cross-sections: " << this->z_step_num << std::endl;
		assert(A.size() == this->z_step_num);

		// change the yarn cross-sections
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];

				const int vertices_num = this->plys[i].fibers[f].vertices.size();
				const int num_of_cores = omp_get_num_procs();
//#pragma omp parallel for num_threads(num_of_cores) 
				for (int v = 0; v < vertices_num; v++) {

					//because of flyaways we should find the closest indx to v
					//int indx = static_cast<int> ((fiber.vertices[v].z - zMin) / this->z_step_size);  //was working before
					int indx = static_cast<int> ( ((fiber.vertices[v].z - zMin) / zSpan) * (this->z_step_num-1)  ); // -1 because index starts from 0

					/* transformation matrix */
					Eigen::Matrix2f transf = A[indx];
					if (global_rot != "") {
						transf = A[indx] * R[indx];
					}

#ifdef SHRINK
					float scale;
					scale = 0.95;
					transf = scale * transf;
#endif

					/* Manually match the phase for short yarns: (6x6 result) */
					//Eigen::Matrix2f globalrot;
					//float theta = -90.0 * pi/180.f;
					//globalrot << cos(theta), -sin(theta), sin(theta), cos(theta);
					//transf = A[indx] * globalrot ;
					/**/


					Eigen::MatrixXf ref(2, 1);
					ref << fiber.vertices[v].x, fiber.vertices[v].y;
					Eigen::MatrixXf def(2, 1);
					def = transf*ref;
					fiber.vertices[v].x = def(0, 0);
					fiber.vertices[v].y = def(1, 0);					
				}
			}
		}
		//plotIntersections("../data/allCrossSection2D_compress.txt", 0.2);
	} // compress_yarn

	void readDeformGradFile(const char* deformGrad, std::vector<Eigen::Matrix3f> &F, std::vector<Eigen::Vector3f> &T) {
		std::ifstream fin;
		fin.open(deformGrad);
		std::string line;

		while (std::getline(fin, line)) {
			std::vector<std::string> splits = split(line, ' ');
			Eigen::Matrix3f deformGrad;
			deformGrad << atof(splits[0].c_str()), atof(splits[1].c_str()), atof(splits[2].c_str()),
				atof(splits[3].c_str()), atof(splits[4].c_str()), atof(splits[5].c_str()),
					atof(splits[6].c_str()), atof(splits[7].c_str()), atof(splits[8].c_str());
			F.push_back(deformGrad);

			//Eigen::Vector3f t;
			//t << atof(splits[9].c_str()), atof(splits[10].c_str()), atof(splits[11].c_str());
			//T.push_back(t);
		}
		fin.close();
	}

	void Yarn::plotIntersections(const char* filename, const float trimPercent = 0.2) {
		// for debuging:
		FILE *fout;
		const int plane_num = this->z_step_num;
		const int f_num = this->plys[0].fibers.size();
		const int ignorPlanes = trimPercent * plane_num; // crop the first and last % of the yarn
		if (fopen_s(&fout, filename, "wt") == 0) {
			fprintf_s(fout, "plane_num: %d \n", plane_num - 2 * ignorPlanes);
			fprintf_s(fout, "ply_num: %d \n \n", this->plys.size());
			for (int step_id = ignorPlanes; step_id < plane_num - ignorPlanes; step_id++) {
				fprintf_s(fout, "index_plane : %d \n", step_id - ignorPlanes);
				for (int i = 0; i < this->plys.size(); i++) {
					fprintf_s(fout, "ply_fiber_num: %d\n", f_num);
					fprintf_s(fout, "plyCenter: %.4f %.4f\n", this->plys[i].fibers[0].vertices[step_id].x, this->plys[i].fibers[0].vertices[step_id].y);
					for (int f = 0; f < f_num; ++f) {
						fprintf_s(fout, "%.4f %.4f\n", this->plys[i].fibers[f].vertices[step_id].x, this->plys[i].fibers[f].vertices[step_id].y);
					}
				}
				fprintf_s(fout, "\n");
			}
			fclose(fout);
		}
	}


	void Yarn::curve_yarn(const char* pntsFile, const char* normsFile, bool scaleXY) {
		std::cout << "step8: map the straight yarn to the spline curve ..." << std::endl;

		/* use hermite spline multiple segments */
		HermiteCurve centerCurve;
		if (normsFile == "")
			centerCurve.init(pntsFile, normsFile);
		else
			//given normals:
			centerCurve.init_norm(pntsFile, normsFile);
		std::vector<Eigen::Vector3d> all_pts, all_tang, all_norm;
		centerCurve.assign(all_pts, all_tang, all_norm);

		double zMin = std::numeric_limits<double>::max(), zMax = std::numeric_limits<double>::lowest();
		for (const auto &ply : plys)
			for (const auto &fiber : ply.fibers)
				for (const auto &vertex : fiber.vertices) {
					zMin = std::min(zMin, static_cast<double>(vertex.z));
					zMax = std::max(zMax, static_cast<double>(vertex.z));
				}
		double zSpan = zMax - zMin;
		double curveLength = centerCurve.totalLength();
		double xyScale = 1.0;

		printf("  zMin: %.4lf, zMax: %.4lf, zSpan: %.4lf\n", zMin, zMax, zSpan);
		printf("  Curve length: %.4lf", curveLength);
		if (scaleXY)
			printf(" (scale: %.4lf)", xyScale = curveLength / zSpan);
		putchar('\n');

		for (auto &ply : plys) {
			for (auto &fiber : ply.fibers) {
				int i = 0;
				for (auto &vertex : fiber.vertices) {

#ifndef IMPROVED_FLYAWAYS
					/* faster but doesn't work for flyaways: */
					Eigen::Vector3d ez = all_tang[i];
					Eigen::Vector3d ey = all_norm[i];
					Eigen::Vector3d ex = ez.cross(ey);

					Eigen::Vector3d pos = all_pts[i];
#else


					double len = curveLength*(vertex.z - zMin) / zSpan;
					double t = centerCurve.arcLengthInvApprox(len);
					// use rotated Frenet frame 
					Eigen::Vector3d ex, ey, ez;
					centerCurve.getRotatedFrame(t, ex, ey, ez);
					Eigen::Vector3d pos = centerCurve.eval(t);
#endif
					Eigen::Vector3d pos1;


					/** local to world **/
					Eigen::Vector3d local;
					local << vertex.x, vertex.y, 0.0; //0 since we are in 2D plane
					Eigen::Matrix3d M;

					M << ex[0], ey[0], ez[0],
						ex[1], ey[1], ez[1],
						ex[2], ey[2], ez[2];

					pos1 = pos + M*local;

					//or:
					//pos1 = pos +xyScale*(static_cast<double>(vertex.x)*ex + static_cast<double>(vertex.y)*ey + 0.0*ez);

					vertex.x = static_cast<float>(pos1[0]);
					vertex.y = static_cast<float>(pos1[1]);
					vertex.z = static_cast<float>(pos1[2]);

					i++;
				}
			}
		}

		//plotIntersections("../data/allCrossSection2D_curve.txt",0.0);
	} // curve_yarn


	void Yarn::write_yarn_obj(const char* filename) {
		printf("Writing vertices ...\n");
		int total_fiber_num = 0, ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++)
			total_fiber_num += this->plys[i].fibers.size();
		std::ofstream fout(filename);
		//fout << total_fiber_num << std::endl; //TODO : generated yarn format should be same as simulated yarn 
		for (int i = 0; i < ply_num; i++) {
			int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				//Fiber &fiber = this->plys[i].fibers[10]; //render one fiber
				Fiber &fiber = this->plys[i].fibers[f];
				int fiber_vertex_num = fiber.vertices.size();
				//fout << fiber_vertex_num << std::endl;
				for (int v = 0; v < fiber_vertex_num; v++) {
					fout << "v " <<  fiber.vertices[v].x << " " << fiber.vertices[v].y << " " << fiber.vertices[v].z << std::endl;
				}
			}
		}
		for (int i = 0; i < ply_num; i++) {
			int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				int fiber_vertex_num = fiber.vertices.size();
				for (int v = 0; v < fiber_vertex_num-1 ; v++) {
					const int indx = i*fiber_num*fiber_vertex_num + f*fiber_vertex_num;
					fout << "l " << indx + v + 1 << " " << indx + v + 2 << std::endl;
				}
			}
		}
		fout.close();
		printf("Writing vertices to OBJ file done!\n");
	}

	void Yarn::write_yarn(const char* filename) {
#ifndef DROPOUT
		printf("Writing vertices ...\n");
		int total_fiber_num = 0, ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++)
			total_fiber_num += this->plys[i].fibers.size();
		std::ofstream fout(filename);
		fout << total_fiber_num << std::endl; //TODO : generated yarn format should be same as simulated yarn 
		for (int i = 0; i < ply_num; i++) {
			int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				//Fiber &fiber = this->plys[i].fibers[10]; //render one fiber
				Fiber &fiber = this->plys[i].fibers[f];
				int fiber_vertex_num = fiber.vertices.size();
				fout << fiber_vertex_num << std::endl;
				for (int v = 0; v < fiber_vertex_num; v++) {
					fout << fiber.vertices[v].x << " " << fiber.vertices[v].y << " " << fiber.vertices[v].z << std::endl;
				}
			}
		}
		fout.close();
		std::cout << "Fibers are written to " << filename << "\n\n";
#else
		const int dropOut = 2; //discard half of the fibers to accelerate rendering process for large textiles
		printf("Writing vertices ...\n");
		int total_fiber_num = 0, ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++)
			total_fiber_num += this->plys[i].fibers.size();
		std::ofstream fout(filename);
		fout << ceil(total_fiber_num / dropOut) << std::endl; //TODO : generated yarn format should be same as simulated yarn 
		for (int i = 0; i < ply_num; i++) {
			int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				if (f % dropOut != 0) continue;
				//Fiber &fiber = this->plys[i].fibers[10]; //render one fiber
				Fiber &fiber = this->plys[i].fibers[f];
				int fiber_vertex_num = fiber.vertices.size();
				fout << fiber_vertex_num << std::endl;
				for (int v = 0; v < fiber_vertex_num; v++) {
					fout << fiber.vertices[v].x << " " << fiber.vertices[v].y << " " << fiber.vertices[v].z << std::endl;
				}
			}
		}
		fout.close();
		printf("Writing vertices to file done!\n");
		std::cout << "\n\n";
#endif

	}


	void Yarn::write_plys(const char *fn) {
#ifdef VERBOSE
		printf("Writing vertices to file...\n");
#endif
		std::string filename;
		int ply_num = this->plys.size();
		//std::string path = WORK_PATH;
		std::vector<std::ofstream> fouts(ply_num);
		for (int i = 0; i < ply_num; i++) {
			if (fn) {
				filename = std::string(fn);
				filename = filename.substr(0, filename.find('.')) + std::to_string(static_cast<long long>(i)) + ".txt";
			}
			else
				filename = (this->output_file.substr(0, this->output_file.find('.')) + std::to_string(static_cast<long long>(i)) + ".txt");
			//filename = (path + this->output_file.substr(0, this->output_file.find('.')) + std::to_string(static_cast<long long>(i)) + ".txt");

			fouts[i].open(filename.c_str());
#ifdef VERBOSE
			std::cout << "Writing File : " << filename << std::endl;
#endif
			int fiber_num = this->plys[i].fibers.size();
#ifndef IMPROVED_FLYAWAYS
			int fly_fiber_num = this->plys[i].fly_fiber_num;
#endif

#ifdef VERBOSE
			std::cout << "fiber_num = " << fiber_num
#ifndef IMPROVED_FLYAWAYS
				<< " flyaway_num = " << fly_fiber_num
#endif
				<< std::endl;
#endif

			fouts[i] << fiber_num
#ifndef IMPROVED_FLYAWAYS
				+ fly_fiber_num
#endif
				<< std::endl;

			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];

				int fiber_vertex_num = fiber.vertices.size();
				fouts[i] << fiber_vertex_num << std::endl;
				for (int v = 0; v < fiber_vertex_num; v++) {
					fouts[i] << fiber.vertices[v].x << " " << fiber.vertices[v].y << " " << fiber.vertices[v].z << std::endl;
				}

#ifndef IMPROVED_FLYAWAYS
				int fly_fiber_num = fiber.fly_vertices_list.size();
				if (fly_fiber_num > 0) {
					for (int fi = 0; fi < fly_fiber_num; fi++) {
						int fly_fiber_vertex = fiber.fly_vertices_list[fi].size();
						fouts[i] << fly_fiber_vertex << std::endl;
						for (int v = 0; v < fly_fiber_vertex; v++)
							fouts[i] << fiber.fly_vertices_list[fi][v].x << " " << fiber.fly_vertices_list[fi][v].y << " " << fiber.fly_vertices_list[fi][v].z << std::endl;
					}
				}
#endif
			}
			fouts[i].close();
		}

#ifdef VERBOSE
		printf("Writing vertices to file done!\n");
#endif
	}


	void Yarn:: L2norm_3D(const Yarn &yarn1, const Yarn &yarn2, const int trimPercent, float &L2) {
		assert(yarn1.getPlyNum() == yarn2.getPlyNum());
		assert(yarn1.getStepNum() == yarn2.getStepNum());
		const int ply_num = yarn1.getPlyNum();

		float e = 0.f;
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = yarn1.plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				const int vertices_num = yarn1.plys[i].fibers[f].vertices.size();
				const int ignorPlanes = trimPercent * vertices_num;
				for (int v = 0; v < vertices_num; v++) {
					const vec3f v1 = yarn1.plys[i].fibers[f].vertices[v];
					const vec3f v2 = yarn2.plys[i].fibers[f].vertices[v];
					if ( v > ignorPlanes && v < vertices_num - ignorPlanes)
						e += square_norm(v1 - v2);
				}
			}
		}
		L2 = e;
	}


} // namespace Fiber