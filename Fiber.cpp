#include "Fiber.h"
#include "hermiteCurve.h"

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
					this->plys[i].fibers.resize(fiber_num);
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


	void Yarn::yarn_simulate() {

		std::cout << "step1: yarn center ...\n";
		const float base_radius = this->yarn_radius;
		this->aabb_procedural.reset();

		std::cout << "step2: initial plys centers ... \n";
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			float angle = 2 * pi * i / ply_num;
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


	void Yarn::compress_yarn(const char* filename) {
		std::cout << "step7: change yarn cross-sections ..." << std::endl;
		// read the compression parameters file
		std::ifstream fin;
		if (filename != NULL)
			fin.open(filename);

		std::string line;

		std::vector<compress> compress_params;
		while (std::getline(fin, line)) {
			compress param;
			std::vector<std::string> splits = split(line, ' ');

			param.z = atof( splits[0].c_str() );
			param.theta = atof( splits[1].c_str() );
			param.ellipse_long = atof( splits[2].c_str() );
			param.ellipse_short = atof( splits[3].c_str() );
			compress_params.push_back(param);
		}
		
		// change the yarn cross-sections
		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				const int vertices_num = this->plys[i].fibers[f].vertices.size();
				for (int v = 0; v < vertices_num; v++) {
					int c = 0; // for interpolation: the z value for the current vertex lies between c-1 and c 
					for (c = 0; c < compress_params.size(); c++) {
						if (v*this->z_step_size < compress_params[c].z) break;
					}

					// linear interpolation to obtain compress-parameters for this z value:
					float ratio = (compress_params[c].theta - compress_params[c - 1].theta) / (compress_params[c].z - compress_params[c - 1].z);
					float theta = ratio * (v*this->z_step_size) + compress_params[c - 1].theta;

					ratio = (compress_params[c].ellipse_long - compress_params[c - 1].ellipse_long) / (compress_params[c].z - compress_params[c - 1].z);
					float ellipse_long = ratio * (v*this->z_step_size) + compress_params[c - 1].ellipse_long;

					ratio = (compress_params[c].ellipse_short - compress_params[c - 1].ellipse_short) / (compress_params[c].z - compress_params[c - 1].z);
					float ellipse_short = ratio * (v*this->z_step_size) + compress_params[c - 1].ellipse_short;

					vec3f axis_x(1.0f, 0.f, 0.f);   //because cross sections in world coord. are defined in xy plane. 
					const float short_axis_x = axis_x.x * std::cosf(theta) - axis_x.y * std::sinf(theta);
					const float short_axis_y = axis_x.x * std::cosf(theta) + axis_x.y * std::sinf(theta);
					vec3f short_axis = vec3f(short_axis_x, short_axis_y, 0.f);
					vec3f long_axis = vec3f(-short_axis.y, short_axis.x, 0);
					float _p_x = nv::dot(fiber.vertices[v], short_axis), _p_y = nv::dot(fiber.vertices[v], long_axis);

					// scale the shape of cross section
					_p_x *= ellipse_short;
					_p_y *= ellipse_long;
					vec3f new_p = _p_x * short_axis + _p_y * long_axis;

					fiber.vertices[v].x = new_p.x;
					fiber.vertices[v].y = new_p.y;
				}
			}
		}
	} // compress_yarn


	void Yarn::curve_yarn(const char* filename, bool scaleXY) {
		std::cout << "step8: map the straight yarn to the spline curve ..." << std::endl;
	
		/* use hermite spline multiple segments */
        HermiteCurve curve;
        curve.init(filename);

        double zMin = std::numeric_limits<double>::max(), zMax = std::numeric_limits<double>::lowest();
        for ( const auto &ply : plys )
            for ( const auto &fiber : ply.fibers )
                for ( const auto &vertex : fiber.vertices ) {
                    zMin = std::min(zMin, static_cast<double>(vertex.z));
                    zMax = std::max(zMax, static_cast<double>(vertex.z));
                }
        double zSpan = zMax - zMin;
        double curveLength = curve.totalLength();
        double xyScale = 1.0;

        printf("  zMin: %.4lf, zMax: %.4lf, zSpan: %.4lf\n", zMin, zMax, zSpan);
        printf("  Curve length: %.4lf", curveLength);
        if ( scaleXY )
            printf(" (scale: %.4lf)", xyScale = curveLength/zSpan);
        putchar('\n');

        for ( auto &ply : plys )
            for ( auto &fiber : ply.fibers )
                for ( auto &vertex : fiber.vertices ) {
                    double len = curveLength*(vertex.z - zMin)/zSpan;
                    double t = curve.arcLengthInvApprox(len);

                    Eigen::Vector3d pos = curve.eval(t), tang = curve.evalTangent(t), norm = curve.evalNormal(t);
                    Eigen::Vector3d binorm = tang.cross(norm);

                    Eigen::Vector3d pos1;
                    pos1 = pos + xyScale*(static_cast<double>(vertex.x)*norm + static_cast<double>(vertex.y)*binorm);

                    vertex.x = static_cast<float>(pos1[0]);
                    vertex.y = static_cast<float>(pos1[1]);
                    vertex.z = static_cast<float>(pos1[2]);
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
	}

} // namespace Fiber