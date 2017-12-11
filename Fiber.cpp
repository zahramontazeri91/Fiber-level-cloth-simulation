#include "Fiber.h"
#include "hermiteCurve.h"
#include "crossSection.h"

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
		//TODO: only initializing "plys", other parameters might be needed to setup later
		std::ifstream fin;
		if (yarnfile != NULL)
			fin.open(yarnfile);

		const int num_of_cores = omp_get_num_procs();

		std::string line;
		std::getline(fin, line);
		int fiber_num = atoi(line.c_str()) / ply_num;
		this->plys.resize(ply_num);
#pragma omp parallel for num_threads(num_of_cores) 
		for (int p = 0; p < ply_num; ++p) {
			this->plys[p].fibers.resize(fiber_num);
			for (int f = 0; f < fiber_num; ++f) {
				Fiber &fiber = this->plys[p].fibers[f];
				fiber.clear(); //clear the vertices list 
				std::getline(fin, line);
				int vrtx_num = atoi(line.c_str());
				for (int v = 0; v < vrtx_num; ++v) {
					std::getline(fin, line);
					std::vector<std::string> splits = split(line, ' ');
					vec3f vrtx(atof(splits[0].c_str()), atof(splits[1].c_str()), atof(splits[2].c_str()));
					fiber.vertices.push_back(vrtx);
				}
			}
		}
		printf("Yarn is initialized from the file. \n");
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

	#pragma omp parallel for num_threads(num_of_cores) 
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
		std::ofstream fout(yarnCenterfile);
		fout << vrtx_num << '\n';
		for (int v = 0; v < vrtx_num; ++v) {
			vec3f sum_fibers = vec3f(0.f);
			for (int f = 0; f < fiber_num; ++f) {
				sum_fibers += yarn.plys[0].fibers[f].vertices[v];
			}
			sum_fibers /= static_cast<float>(fiber_num);
			fout << sum_fibers.x << " " << sum_fibers.y << " " << sum_fibers.z << "\n";
		}
		fout.close();
		
		printf("Centerfiber is written to the file. \n");
	}

	void Yarn::parse(const char* filename) {
		std::ifstream fin;
		if (filename != NULL)
			fin.open(filename);

		std::string line;

		while (std::getline(fin, line)) {
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
					this->plys[i].fibers.resize(fiber_num + 1); //add one fiber for ply-center 
				}
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
		std::cout << "Parsing the input is done!\n\n";
	}


	Yarn::Yarn() {}
	Yarn::~Yarn() {}

	void Yarn::assignParameterizePlyCenters(const char *plyCenterFile) {
		std::ifstream fin;
		fin.open(plyCenterFile);
		const int ply_num = this->plys.size();
		std::string line;
		std::getline(fin, line);
		float plane_num = atoi(line.c_str());
		/* initialize first fiber for each ply as the ply-center */
		for (int step_id = 0; step_id < plane_num; step_id++) {
			const float z = this->z_step_size * (step_id - this->z_step_num / 2.f); // devided by 2 Bcuz yarn lies between neg and pos z
			std::string line;
			std::getline(fin, line);
			std::vector<std::string> splits = split(line, ' ');

			float R = atof(splits[0].c_str());
			float theta = atof(splits[1].c_str());

			for (int i = 0; i < ply_num; i++) {

				float world_x = R * cos(theta + i * 2 * pi / ply_num);
				float world_y = R * sin(theta + i * 2 * pi / ply_num);

				vec3f verIn = vec3f(world_x, world_y, z), verOut;
				verOut = verIn;

				this->aabb_procedural.grow(verOut);
				if (this->aabb_micro_ct.in(verOut))
					this->plys[i].fibers[0].vertices.push_back(verOut);
			}
		}
		fin.close();
	}

	void Yarn::assignPlyCenters(const char *plyCenterFile) {
		std::ifstream fin;
		fin.open(plyCenterFile);
		const int ply_num = this->plys.size();

		/* initialize first fiber for each ply as the ply-center */
		for (int step_id = 0; step_id < this->z_step_num; step_id++) {
			for (int i = 0; i < ply_num; i++) {
				const float z = this->z_step_size * (step_id - this->z_step_num / 2.f); // devided by 2 Bcuz yarn lies between neg and pos z
				std::string line;
				std::getline(fin, line);
				std::vector<std::string> splits = split(line, ' ');
				float world_x = atof(splits[0].c_str());
				float world_y = atof(splits[1].c_str());

				vec3f verIn = vec3f(world_x, world_y, z), verOut;
				verOut = verIn;

				this->aabb_procedural.grow(verOut);
				if (this->aabb_micro_ct.in(verOut))
					this->plys[i].fibers[0].vertices.push_back(verOut);
			}
			std::string whitespace;
			std::getline(fin, whitespace);
		}
		fin.close();
	}

	float Yarn::extractFiberTwist(const char *fiberTwistFile) {
		std::ifstream fin;
		if (fiberTwistFile != NULL)
			fin.open(fiberTwistFile);

		float fiberTwistAVG = 0.f;
		std::string line;
		std::getline(fin, line);
		int plane_num = atoi(line.c_str());
		const int ignorPlanes = 0.1 * plane_num; // crop the first and last 10% of the yarn
		int i = 0;
		for (i = 0; i < plane_num; ++i) {
			std::getline(fin, line);
			if (i > ignorPlanes && i < plane_num - ignorPlanes) {
				float fiberTwist = atof(line.c_str());
				fiberTwistAVG += fiberTwist;
			}
		}
		fiberTwistAVG /= static_cast<float>(plane_num - 2 * ignorPlanes);
		return fiberTwistAVG;
	}

	void Yarn::yarn_simulate(const char *plyCenterFile, const char *fiberTwistFile) {

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

		std::cout << "step4-5-6: rotate ply-centers around yarn-center and fibers around ply-centers and apply the compression ... \n";
#pragma omp parallel for num_threads(num_of_cores)
		/* initialize the yarn fibers by assigning the ply-center to the fiber[0]*/
		//assignPlyCenters(plyCenterFile);
		assignParameterizePlyCenters(plyCenterFile);
		/* Calculate the fiber twisting by averaging each fiber rotation in two consecutive cross-section and average over all planes */
		const float fiber_theta_avg = extractFiberTwist(fiberTwistFile);

		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
			// generate all fibers around ply-center
			for (int f = 1; f < fiber_num; f++) { //starts from index 1 because index0 is reserved for ply-center
				Fiber &fiber = this->plys[i].fibers[f];
				fiber.clear(); //clear the vertices list 
				for (int step_id = 0; step_id < this->z_step_num; step_id++) {
					const float z = this->z_step_size * (step_id - this->z_step_num / 2.f); // devided by 2 Bcuz yarn lies between neg and pos z
																							//const float fiber_theta = this->plys[i].clock_wise ? -z * 2 * pi / this->plys[i].alpha : z * 2 * pi / this->plys[i].alpha;
					const float fiber_theta = fiber_theta_avg;

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

					// translate it to ply-center (fiber_0)
					float plyCenter_x = this->plys[i].fibers[0].vertices[step_id].x;
					float plyCenter_y = this->plys[i].fibers[0].vertices[step_id].y;
					world_x = local_x + plyCenter_x;
					world_y = local_y + plyCenter_y;

					//render only the plyCenters
					//world_x = plyCenter_x;
					//world_y = plyCenter_y;

					vec3f verIn = vec3f(world_x, world_y, z), verOut;
					verOut = verIn;

					this->aabb_procedural.grow(verOut);
					if (this->aabb_micro_ct.in(verOut))
						fiber.vertices.push_back(verOut);
				}
			}
		}
		/*
		// for debuging:
		FILE *fout;
		if (fopen_s(&fout, "../data/plyCenter_proc.txt", "wt") == 0) {
		for (int step_id = 0; step_id < this->z_step_num; step_id++) {
		for (int i = 0; i < ply_num; i++) {
		fprintf_s(fout, "%.6f %.6f\n", this->plys[i].fibers[0].vertices[step_id].x, this->plys[i].fibers[0].vertices[step_id].y);
		}
		fprintf_s(fout, "\n");
		}
		fclose(fout);
		}
		*/
		FILE *fout;
		// for debuging:
		const int plane_num = this->z_step_num;
		const int f_num = this->plys[0].fibers.size() - 1; //since the first fiber is the center
		const int ignorPlanes = 0.1 * plane_num; // crop the first and last 10% of the yarn
		if (fopen_s(&fout, "../data/allCrossSection2D_proc.txt", "wt") == 0) {
			fprintf_s(fout, "plane_num: %d \n", plane_num - 2 * ignorPlanes);
			fprintf_s(fout, "ply_num: %d \n \n", ply_num);
			for (int step_id = ignorPlanes; step_id < plane_num - ignorPlanes; step_id++) {
				for (int i = 0; i < ply_num; i++) {
					fprintf_s(fout, "ply_fiber_num: %d\n", f_num);
					fprintf_s(fout, "plyCenter: %.4f %.4f\n", this->plys[i].fibers[0].vertices[step_id].x, this->plys[i].fibers[0].vertices[step_id].y);
					for (int f = 1; f < f_num + 1; ++f) {
						fprintf_s(fout, "%.4f %.4f\n", this->plys[i].fibers[f].vertices[step_id].x, this->plys[i].fibers[f].vertices[step_id].y);
					}
				}
				fprintf_s(fout, "\n");
			}
			fclose(fout);
		}
	} // yarn_simulate

#if 0
	void Yarn::yarn_simulate() {

		std::cout << "step1: yarn center ...\n";
		const float base_radius = this->yarn_radius;
		this->aabb_procedural.reset();

		std::cout << "step2: initial plys centers ... \n";
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num; // place ply centers in a circle with radius yarn-radius/2 

												//TODO: this works well for 2-ply yarn only. How about 5-ply yarn?
			this->plys[i].base_theta = angle;
			this->plys[i].base_center = vec3f(base_radius / 2 * std::cosf(angle), base_radius / 2 * std::sinf(angle), 0);
		}

		const int num_of_cores = omp_get_num_procs();

		std::cout << "step3: initial plys locations ... \n";
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

		std::cout << "step4-5-6: rotate ply-centers around yarn-center and fibers around ply-centers and apply the compression ... \n";
#pragma omp parallel for num_threads(num_of_cores)

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
#endif

	void Yarn::readCompressFile(const char* filename, std::vector<compress> &compress_params) {
		std::ifstream fin;
		if (filename != NULL)
			fin.open(filename);

		std::string line;
		std::getline(fin, line);
		const int plane_num = atof(line.c_str());
		int planeId = 0;
		while (std::getline(fin, line)) {
			compress param;
			std::vector<std::string> splits = split(line, ' ');

			param.ellipse_long = atof(splits[0].c_str());
			param.ellipse_short = atof(splits[1].c_str());
			param.theta = atof(splits[2].c_str());
			compress_params.push_back(param);
			planeId++;
		}
	}


	void Yarn::fitCircle(const yarnIntersect2D &pts, float &radius)
	{
		//Find the total number of points for all plys
		int sz = 0;
		for (int p = 0; p < pts.size(); ++p)
			sz += pts[p].size();
		assert(sz);
		cv::Mat data_pts(sz, 2, CV_32FC1, cv::Scalar::all(0));

		int c = data_pts.rows;
		for (int p = 0; p < pts.size(); ++p) {
			for (int i = 0; i < pts[p].size(); ++i) {

				data_pts.at<float>(c, 0) = pts[p][i].x;
				data_pts.at<float>(c, 1) = pts[p][i].y;
				--c;
			}
		}

		//Perform PCA analysis
		cv::PCA pca_analysis(data_pts, cv::Mat(), cv::PCA::DATA_AS_ROW, 2);

		//Store the center of the object
		vec2f center = vec2f(pca_analysis.mean.at<float>(0, 0), pca_analysis.mean.at<float>(0, 1));
		//Store the eigenvalues and eigenvectors
		std::vector<vec2f> eigen_vecs(1);
		//std::vector<float> eigen_val(1);
		eigen_vecs[0] = vec2f(pca_analysis.eigenvectors.at<float>(0, 0),
			pca_analysis.eigenvectors.at<float>(0, 1));

		//find eigen values by projecting the points on eigenVectors
		float max = std::numeric_limits<float>::min();
		for (int i = 0; i < data_pts.rows; ++i) {
			float prj = eigen_vecs[0].x * data_pts.at<float>(i, 0) + eigen_vecs[0].y * data_pts.at<float>(i, 1);
			if (std::abs(prj) > max)
				max = std::abs(prj);
		}
		radius = max;
	}


	void Yarn::yarn2crossSections(std::vector<yarnIntersect2D> &itsLists) {
		//first initialize the vectors
		itsLists.resize(this->z_step_num);
		for (int i = 0; i < itsLists.size(); ++i)
			itsLists[i].resize(this->plys.size());

		//copy the yarn into new dataStructure
		for (int p = 0; p < this->plys.size(); ++p) {
			plyItersect plyIts;
			for (int f = 0; f < this->plys[p].fibers.size(); ++f) {
				for (int v = 0; v < this->plys[p].fibers[f].vertices.size(); ++v) {
					itsLists[v][p].push_back(vec2f(this->plys[p].fibers[f].vertices[v].x,
						this->plys[p].fibers[f].vertices[v].y));
				}
			}
		}
	}

	void Yarn::compress_yarn(const char* filename) {
		std::cout << "step7: compress yarn cross-sections ..." << std::endl;

		//first find the fitted circle around each cross-section
		std::vector<float> fitCircleR;
		const int ply_num = this->plys.size();
		std::vector<yarnIntersect2D> itsLists;
		yarn2crossSections(itsLists);
		float fitCircleR_avg = 0.f;
		for (int i = 0; i < this->z_step_num; ++i) {
			float radius;
			fitCircle(itsLists[i], radius);
			fitCircleR.push_back(radius);
			fitCircleR_avg += radius;
		}
		fitCircleR_avg /= static_cast<float> (this->z_step_num);

		std::vector<compress> compress_params;
		readCompressFile(filename, compress_params);

		// change the yarn cross-sections
		//const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				const int vertices_num = this->plys[i].fibers[f].vertices.size();
				//assert(compress_params.size() == vertices_num);
				//std::cout << compress_params.size() << " " <<  vertices_num << std::endl;
				//while (1);
				for (int v = 0; v < vertices_num; v++) {

					const float ellipse_long = compress_params[v].ellipse_long;
					const float ellipse_short = compress_params[v].ellipse_short;
					const float ellipse_theta = compress_params[v].theta;

					//obtain ellipse axis
					const vec2f ellipse_axis_long = vec2f(cos(ellipse_theta), sin(ellipse_theta));
					const vec2f ellipse_axis_short = vec2f(-sin(ellipse_theta), cos(ellipse_theta));

					//transfer from x-y space to ellipse space
					vec2f world_p(fiber.vertices[v].x, fiber.vertices[v].y);
					vec2f ellipse_p(0.f);
					ellipse_p.x = nv::dot(ellipse_axis_long, world_p);
					ellipse_p.y = nv::dot(ellipse_axis_short, world_p);

					//apply the scaling 
					ellipse_p.x *= ellipse_long / fitCircleR_avg;
					ellipse_p.y *= ellipse_short / fitCircleR_avg;
					//ellipse_p.x *= ellipse_long / fitCircleR[v];
					//ellipse_p.y *= ellipse_short / fitCircleR[v];

					//transfer back to x-y
					world_p.x = nv::dot(vec2f(ellipse_axis_long.x, ellipse_axis_short.x), ellipse_p);
					world_p.y = nv::dot(vec2f(ellipse_axis_long.y, ellipse_axis_short.y), ellipse_p);

					fiber.vertices[v].x = world_p.x;
					fiber.vertices[v].y = world_p.y;
				}
			}
		}
		// for debuging:
		FILE *fout;
		const int plane_num = this->z_step_num;
		const int f_num = this->plys[0].fibers.size() - 1; //since the first fiber is the center
		const int ignorPlanes = 0.1 * plane_num; // crop the first and last 10% of the yarn
		if (fopen_s(&fout, "../data/allCrossSection2D_compress.txt", "wt") == 0) {
			fprintf_s(fout, "plane_num: %d \n", plane_num - 2 * ignorPlanes);
			fprintf_s(fout, "ply_num: %d \n \n", ply_num);
			for (int step_id = ignorPlanes; step_id < plane_num - ignorPlanes; step_id++) {
				for (int i = 0; i < ply_num; i++) {
					fprintf_s(fout, "ply_fiber_num: %d\n", f_num);
					fprintf_s(fout, "plyCenter: %.4f %.4f\n", this->plys[i].fibers[0].vertices[step_id].x, this->plys[i].fibers[0].vertices[step_id].y);
					for (int f = 1; f < f_num + 1; ++f) {
						fprintf_s(fout, "%.4f %.4f\n", this->plys[i].fibers[f].vertices[step_id].x, this->plys[i].fibers[f].vertices[step_id].y);
					}
				}
				fprintf_s(fout, "\n");
			}
			fclose(fout);
		}
	} // compress_yarn


	void Yarn::curve_yarn(const char* pntsFile, const char* normsFile, bool scaleXY) {
		std::cout << "step8: map the straight yarn to the spline curve ..." << std::endl;

		/* use hermite spline multiple segments */
		HermiteCurve curve;
		//curve.init(pntsFile);
		curve.init(pntsFile, normsFile);

		double zMin = std::numeric_limits<double>::max(), zMax = std::numeric_limits<double>::lowest();
		for (const auto &ply : plys)
			for (const auto &fiber : ply.fibers)
				for (const auto &vertex : fiber.vertices) {
					zMin = std::min(zMin, static_cast<double>(vertex.z));
					zMax = std::max(zMax, static_cast<double>(vertex.z));
				}
		double zSpan = zMax - zMin;
		double curveLength = curve.totalLength();
		double xyScale = 1.0;

		printf("  zMin: %.4lf, zMax: %.4lf, zSpan: %.4lf\n", zMin, zMax, zSpan);
		printf("  Curve length: %.4lf", curveLength);
		if (scaleXY)
			printf(" (scale: %.4lf)", xyScale = curveLength / zSpan);
		putchar('\n');

		//
		//bool firstTime = true;
		//std::ofstream fout1("genYarn_frame29_norms.txt");
		//fout1 << this->z_step_num << "\n";
		//

		for (auto &ply : plys)
			for (auto &fiber : ply.fibers) {
				for (auto &vertex : fiber.vertices) {
					double len = curveLength*(vertex.z - zMin) / zSpan;
					double t = curve.arcLengthInvApprox(len);

					// use rotated Frenet frame 
					Eigen::Vector3d ex, ey, ez;
					curve.getRotatedFrame(t, ex, ey, ez);

					Eigen::Vector3d pos = curve.eval(t);
					Eigen::Vector3d pos1;
					pos1 = pos + xyScale*(static_cast<double>(vertex.x)*ex + static_cast<double>(vertex.y)*ey);

					vertex.x = static_cast<float>(pos1[0]);
					vertex.y = static_cast<float>(pos1[1]);
					vertex.z = static_cast<float>(pos1[2]);

					//for debugging (write normals to file):
					//if (firstTime) {
					//fout1 << ey[0] << " " << ey[1] << " " << ey[2] << std::endl;
					//}
				}
		//firstTime = false;
		//fout1.close();
		}

		// for debuging:
		const int ply_num = this->plys.size();
		FILE *fout;
		const int plane_num = this->z_step_num;
		const int f_num = this->plys[0].fibers.size() - 1; //since the first fiber is the center
		const int ignorPlanes = 0.1 * plane_num; // crop the first and last 10% of the yarn
		if (fopen_s(&fout, "../data/allCrossSection2D_curve.txt", "wt") == 0) {
			fprintf_s(fout, "plane_num: %d \n", plane_num - 2 * ignorPlanes);
			fprintf_s(fout, "ply_num: %d \n \n", ply_num);
			for (int step_id = ignorPlanes; step_id < plane_num - ignorPlanes; step_id++) {
				for (int i = 0; i < ply_num; i++) {
					fprintf_s(fout, "ply_fiber_num: %d\n", f_num);
					fprintf_s(fout, "plyCenter: %.4f %.4f\n", this->plys[i].fibers[0].vertices[step_id].x, this->plys[i].fibers[0].vertices[step_id].y);
					for (int f = 1; f < f_num + 1; ++f) {
						fprintf_s(fout, "%.4f %.4f\n", this->plys[i].fibers[f].vertices[step_id].x, this->plys[i].fibers[f].vertices[step_id].y);
					}
				}
				fprintf_s(fout, "\n");
			}
			fclose(fout);
		}
	} // curve_yarn


	void Yarn::write_yarn(const char* filename) {
		std::cout << "\n\n";
		printf("Writing vertices ...\n");
		int total_fiber_num = 0, ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++)
			total_fiber_num += this->plys[i].fibers.size();
		std::ofstream fout(filename);
		fout << total_fiber_num << std::endl; //TODO : generated yarn format should be same as simulated yarn 
		for (int i = 0; i < ply_num; i++) {
			int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
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

		//for debugging:
		std::ofstream fout0("genYarn_ply0.txt");
		std::ofstream fout1("genYarn_ply1.txt");
		fout0 << total_fiber_num << std::endl; //TODO : generated yarn format should be same as simulated yarn 
		fout1 << total_fiber_num << std::endl;
		for (int i = 0; i < ply_num; i++) {
			int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				int fiber_vertex_num = fiber.vertices.size();
				if (i == 0)
					fout0 << fiber_vertex_num << std::endl;
				else
					fout1 << fiber_vertex_num << std::endl;
				for (int v = 0; v < fiber_vertex_num; v++) {
					if (i == 0)
						fout0 << fiber.vertices[v].x << " " << fiber.vertices[v].y << " " << fiber.vertices[v].z << std::endl;
					else
						fout1 << fiber.vertices[v].x << " " << fiber.vertices[v].y << " " << fiber.vertices[v].z << std::endl;
				}
			}
		}
		fout0.close();
		fout1.close();

	}
	//
	//void Yarn::shapeCrossSection(yarnIntersect2D &its, float &rLong, float &rShort) {

	//	//Find the total number of points for all plys
	//	int sz = 0;
	//	std::vector<vec2f> pnts;
	//	for (int p = 0; p < its.size(); ++p) {
	//		sz += its[p].size();
	//		for (int v = 0; v < its[p].size(); ++v) {
	//			pnts.push_back( its[p][v] );
	//		}
	//	}


	//	//find eigen values by projecting the points on eigenVectors
	//	float max0 = std::numeric_limits<float>::min();
	//	float max1 = std::numeric_limits<float>::min();
	//	for (int i = 0; i < sz; ++i) {
	//		if (std::abs(pnts[i].x) > max0)
	//			max0 = std::abs(pnts[i].x);
	//		if (std::abs(pnts[i].y) > max1)
	//			max1 = std::abs(pnts[i].y);
	//	}
	//	rLong = std::max(max0,max1);
	//	rShort = std::min(max0,max1);
	//}

} // namespace Fiber