#include "Fiber.h"
#include "hermiteSeg.h"
//#include <Eigen/Dense>

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
					//// std::cout << compress_params[c-1].z << "  " << v*this->z_step_size << "  " << compress_params[c].z << std::endl;
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
	}

	void Yarn::curve_yarn(const char* filename) {
		std::cout << "step8: map the straight yarn to the spline curve ..." << std::endl;
	
		/* use hermite spline multiple segments */
		HermiteSpline_multiSeg splines(filename);

		/* Hermite spline for one-segment */
		//const double first_z = this->plys[0].fibers[0].vertices[0].z;
		//const double last_z = this->plys[0].fibers[0].vertices[this->plys[0].fibers[0].vertices.size() - 1].z;
		//Eigen::Vector3d start_p(0, 0, first_z), end_p(0, 0, last_z);
		//Eigen::Vector3d start_tg(0,1,0), end_tg(0,-1,0);
		//HermiteSpline spline(start_p, end_p, start_tg, end_tg );		

		/* using the paper "Calculation of reference frames along a space curve" */
		vec3f T0, N0, B0;
		vec3f T1, N1, B1;
		Eigen::Vector3d world;

		std::ofstream fout_p ("spline_positions.txt");
		std::ofstream fout_t ("spline_tangents.txt");
		std::ofstream fout_n ("spline_normals.txt");
		bool file_done = false; //write to file only for the first fiber

		const int ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++) {
			const int fiber_num = this->plys[i].fibers.size();
			for (int f = 0; f < fiber_num; f++) {
				Fiber &fiber = this->plys[i].fibers[f];
				const int vertices_num = this->plys[i].fibers[f].vertices.size();
				for (int v = 0; v < vertices_num; v++) {

					/* map to a sin wave (1/10 sin (10x) )  */
					//const float p_z = fiber.vertices[v].z;
					//const vec3f P (0.f, (0.1 * std::sinf(10.0*p_z)), p_z);
					//const vec3f V (0.f, std::cosf(10.0 * p_z), p_z); //derivative of the spline
					//const vec3f Q (0.f, -10.0 * std::sinf(10.0 * p_z), p_z); //derivative of the V (velocity)

					/* map to a given spline */				
					//double t = splines.get_seg_num() * curve / curve_total;
					double t = splines.get_seg_num() * double(v) / double(vertices_num);

					Eigen::Vector3d P_e = splines.eval(t);
					Eigen::Vector3d V_e = splines.evalTangent(t);
					Eigen::Vector3d Q_e = splines.evalCurvature(t);
					const vec3f P = vec3f(P_e[0], P_e[1], P_e[2]);
					const vec3f V = vec3f(V_e[0], V_e[1], V_e[2]);
					const vec3f Q = vec3f(Q_e[0], Q_e[1], Q_e[2]);				

					 
					//if (t > 0.99 && t < 1.01)
						//std::cout << t << " *************\n";
					if ( !(t-int(t)) ) { // TODO: Initialize the Frenet frame for the first point of each spline
						// obtain T, N, and B vectors for the first cross section of the yarn						
						T0 = nv::normalize(V);
						if (Q == vec3f(0.f))
							std::cout << "error: If curvature is zero, N can be any perpendiculat vector to T. \n";
						//N0 = nv::normalize(cross(cross(V, Q), V));
						N0 = vec3f(1, 0, 0);
						B0 = cross(T0, N0);
						
						//use matrix for transformation
						Eigen::Vector3d T_e(T0.x, T0.y, T0.z);
						Eigen::Vector3d N_e(N0.x, N0.y, N0.z);
						Eigen::Vector3d B_e(B0.x, B0.y, B0.z);
						Eigen::MatrixXd R(3, 3);
						R << N_e, B_e, T_e;
						Eigen::Vector3d local(fiber.vertices[v].x, fiber.vertices[v].y, 0.f);
						world = R*local;

						if (!file_done /*&& !(v % 5)*/ ) //write one third 
						{
							fout_p << P.x << " " << P.y << " " << P.z << std::endl;
							fout_t << T0.x << " " << T0.y << " " << T0.z << std::endl;
							fout_n << N0.x << " " << N0.y << " " << N0.z << std::endl;
						}
					}
					/* find N and B for the subsequence cross-sections */
					else {
						T1 = nv::normalize(V);
						vec3f A = cross(T0, T1) / ( nv::length(T0)*nv::length(T1) );
						float sqx = A.x * A.x;
						float sqy = A.y * A.y;
						float sqz = A.z * A.z;
						float cos = dot(T0, T1);
						if (cos >= 1.f) { // If the tangent does not change, nor should the normal
							B1 = B0;
							N1 = N0;
							continue;
						}

						float cos1 = 1.f - cos;
						float xycos1 = A.x * A.y * cos1;
						float yzcos1 = A.y * A.z * cos1;
						float zxcos1 = A.x * A.z * cos1;
						float sin = sqrt(1 - cos*cos);
						float xsin = A.x * sin;
						float ysin = A.y * sin;
						float zsin = A.z * sin;

						Eigen::MatrixXd R_frenet(3, 3);
						R_frenet << sqx + (1 - sqx)*cos,    xycos1 + zsin,           zxcos1 - ysin,
									xycos1 - zsin,           sqy + (1 - sqy)*cos,     yzcos1 + xsin,
									zxcos1 + ysin,           yzcos1 - xsin,           sqz + (1 - sqz)*cos;		

						Eigen::Vector3d B0_e (B0.x, B0.y, B0.z);
						Eigen::Vector3d N0_e (N0.x, N0.y, N0.z);
						Eigen::Vector3d B1_e = R_frenet*B0_e;
						Eigen::Vector3d N1_e = R_frenet*N0_e;
						//Eigen::Vector3d N1_e(1.0, 0.0, 0.0);
						Eigen::Vector3d T1_e (T1.x, T1.y, T1.z);
						B1 = vec3f(B1_e[0], B1_e[1], B1_e[2]);
						N1 = vec3f(N1_e[0], N1_e[1], N1_e[2]);

						//normalize the new vectors
						B1 = nv::normalize(B1);
						N1 = nv::normalize(N1);
						N1_e = Eigen::Vector3d(B1.x, B1.y, B1.z);
						N1_e = Eigen::Vector3d(N1.x, N1.y, N1.z);

						//use matrix for transformation
						Eigen::MatrixXd R(3, 3);
						R << N1_e, B1_e, T1_e;
						Eigen::Vector3d local(fiber.vertices[v].x, fiber.vertices[v].y, 0.f);
						world = R*local;

						T0 = T1;
						B0 = B1;
						N0 = N1;

						// DEBUG: write results to file
						if (!file_done)// && !(v % 5)) //write one third 
						{
							fout_p << P.x << " " << P.y << " " << P.z << std::endl;
							fout_t << T1.x << " " << T1.y << " " << T1.z << std::endl;
							fout_n << N1.x << " " << N1.y << " " << N1.z << std::endl;
						}
					}
					fiber.vertices[v] = P + vec3f(world[0], world[1], world[2]);
				}
				file_done = true;

				//std::cout << "Fiber " << f << " of ply " << i << " is generated.\n";
			}
		}
		fout_p.close();
		fout_t.close();
		fout_n.close();
	}

	void Yarn::write_yarn(const char* filename) {
		std::cout << "\n\n";
		printf("Writing vertices ...\n");
		int total_fiber_num = 0, ply_num = this->plys.size();
		for (int i = 0; i < ply_num; i++)
			total_fiber_num += this->plys[i].fibers.size(); 
		std::ofstream fout(filename);
		fout << total_fiber_num << std::endl;

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