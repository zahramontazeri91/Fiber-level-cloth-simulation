#include "hermiteTests.h"
#include "crossSectionTests.h"
#include "../crossSection.h"

using namespace std;

void linePlaneIntersection_test() {

	const char* yarnfile = "genYarn.txt";
	const char* curvefile = "frame00029_avg.txt";
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);

	Plane plane;
	plane.n = vec3f(1, 0, 0);
	plane.point = vec3f(0, 0, 0);
	vec3f its(0.f);
	vec3f start(1, 2, 0);
	vec3f end(-1, 0, 0);
	bool a = cs.linePlaneIntersection(start, end, plane, its);
	std::cout << " Intersection: " << a << "  " << its.x << " " << its.y << " " << its.z << std::endl;
}

void yarnPlaneIntersection_test() {
	const char* yarnfile = "genYarn_1.txt";
	const char* curvefile = "frame00029_avg.txt"; 
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);

	Plane plane;
	cs.get_plane(1, plane);

	//Plane plane;
	//plane.normal = vec3f(-0.0316426, -0.304221,  0.952076);
	//plane.point = vec3f(0,0, plane0.point.z);

	std::cout << plane.n.x << "  " << plane.n.y << "  " << plane.n.z << std::endl;
	std::cout << plane.point.x << "  " << plane.point.y << "  " << plane.point.z << std::endl;

	yarnIntersect itsList;
	bool b = cs.yarnPlaneIntersection(plane, itsList);

	//std::cout << itsList[0].size() << std::endl;

	plyItersect ply_its = itsList[0];
	//write for only one ply
	FILE *fout;
	if (fopen_s(&fout, "../data/crossSection.txt", "wt") == 0) {
		fprintf_s(fout, "%d \n", ply_its.size());
		// First write the yarn-center
		fprintf_s(fout, "%.4lf %.4lf %.4lf \n", plane.point[0], plane.point[1], plane.point[2]);
		for (int i = 0; i < ply_its.size(); ++i) {
			fprintf_s(fout, "%.4lf %.4lf %.4lf \n", ply_its[i][0], ply_its[i][1], ply_its[i][2]);
		}
		fclose(fout);
	}
}

void buildPlanes_test() {
	const char* yarnfile = "genYarn_1.txt";
	const char* curvefile = "frame00029_avg.txt";
	const int num_planes = 1000;
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);
	FILE *fout;
	if (fopen_s(&fout, "../data/test_planes.txt", "wt") == 0) {
		fprintf_s(fout, "%d \n", num_planes);
		for (int i = 0; i < num_planes; i++) {
			Plane plane;
			cs.get_plane(i, plane);
			fprintf_s(fout, "%.4lf %.4lf %.4lf \n", plane.n.x, plane.n.y, plane.n.z);
		}
		fclose(fout);
	}
}

void allPlanesIntersections_test() {
	//const char* yarnfile = "genYarn.txt"; //For procedural yarn
	const char* yarnfile = "frame00029_scaled.txt"; //For simulated yarn
	const char* curvefile = "frame00029_avg.txt"; 
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	//std::cout << "intersections lists size: " << itsLists.size() << std::endl;

	cs.write_PlanesIntersections3D("../data/allCrossSection.txt",itsLists);
}

void project2Plane_test() {
	const char* yarnfile = "frame00001_scaled.txt"; //For simulated yarn
	const char* curvefile = "frame00029_avg.txt"; 
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);

	vec3f pnt = itsLists[4][0][2];
	vec2f prj;
	Plane plane;
	cs.get_plane(4, plane);
	cs.project2Plane(pnt, plane, prj );
	std::cout << "3D point is: "<< pnt.x << " " << pnt.y << " " << pnt.z << std::endl;
	std::cout << "Projecte point is: " << prj.x << " " << prj.y << std::endl;
}

void write_PlanesIntersections2D_test() {
	const char* yarnfile = "genYarn.txt"; //For procedural yarn
	//const char* yarnfile = "frame00029_scaled.txt"; //For simulated yarn
	const char* curvefile = "frame00029_avg.txt"; 
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);

	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	std::cout << "intersections lists size: " << itsLists.size() << std::endl;

	//std::vector<yarnIntersect2D> deCompressPlaneIntersect;
	//cs.deCompressYarn(allPlaneIntersect, ellipses, deCompressPlaneIntersect);
	//allPlaneIntersect = deCompressPlaneIntersect;

	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);
	std::cout << "number of intersected planes: " << allPlaneIntersect.size() << std::endl;

	//write the 2D intersections for tetsting
	FILE *fout;
	if (fopen_s(&fout, "../data/allCrossSection2D.txt", "wt") == 0) {

		const int ignorPlanes = 0.1 * allPlaneIntersect.size(); // crop the first and last 10% of the yarn
		fprintf_s(fout, "plane_num: %d \n", allPlaneIntersect.size() - 2 * ignorPlanes);
		fprintf_s(fout, "ply_num: %d \n", allPlaneIntersect[0].size());
		fprintf_s(fout, "\n");

		for (int i = ignorPlanes; i < allPlaneIntersect.size() - ignorPlanes; ++i) { //number of planes

			//Plane plane;
			//cs.get_plane(i, plane);
			//fprintf_s(fout, "center: %.4lf %.4lf %.4lf \n", plane.point[0], plane.point[1], plane.point[2]); //center is always 0,0
			for (int p = 0; p < allPlaneIntersect[i].size(); ++p) { //number of plys
				fprintf_s(fout, "ply_fiber_num: %d \n", allPlaneIntersect[i][p].size());

				//find ply center
				vec2f plyCenter(0.f);
				for (int j = 0; j < allPlaneIntersect[i][p].size(); ++j) { //number of intersections
					plyCenter += allPlaneIntersect[i][p][j];
				}
				plyCenter /= allPlaneIntersect[i][p].size();
				fprintf_s(fout, "plyCenter: %.4lf %.4lf \n", plyCenter.x, plyCenter.y);
				for (int j = 0; j < allPlaneIntersect[i][p].size(); ++j) { //number of intersections
					fprintf_s(fout, "%.4f %.4f \n", allPlaneIntersect[i][p][j].x, allPlaneIntersect[i][p][j].y);
				}
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
}

void getOrientation_test() {
	//const char* yarnfile = "gen_yarn_f1.txt"; //For procedural yarn
	const char* yarnfile = "frame00029_scaled.txt"; //For simulated yarn
	const char* curvefile = "frame00029_avg.txt";
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);

	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);

	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);

	//transfer from e1-e2 to x-y plane
	std::vector<yarnIntersect2D> xy_Its;
	cs.transferLocal2XY(allPlaneIntersect, xy_Its);

	Ellipse e;
	cs.fitEllipse(xy_Its[28], e);
	FILE *fout;
	if (fopen_s(&fout, "../data/pca_test.txt", "wt") == 0) {
		for (int p = 0; p < xy_Its[28].size(); ++p) {
			fprintf_s(fout, "%d \n", xy_Its[28][p].size());
			for (int i = 0; i < xy_Its[28][p].size(); ++i) {
				fprintf_s(fout, "%.6f %.6f \n", xy_Its[28][p][i].x, xy_Its[28][p][i].y);
			}
		}
		fprintf_s(fout, "\n");
		fprintf_s(fout, "%.4f %.4f \n", e.center.x, e.center.y);
		fprintf_s(fout, "%.4f %.4f %.4f \n", e.longR, e.shortR, e.angle);
		fprintf_s(fout, "\n");
		fclose(fout);
	}
}
void extractCompressParam_test() {
	/* ply-center simulate: 
	 write avg of each ply */
	/* ply-center procedural:
	read generated-yarn and write fiber[0] */

	const char* yarnfile = "genYarn.txt"; //For procedural yarn
	//const char* yarnfile = "frame00029_scaled.txt"; //For simulated yarn
	const char* curvefile = "frame00029_avg.txt";
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);

	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect, ellipses);

	// Decompress simulated yarn
	//std::vector<yarnIntersect2D> deCompressPlaneIntersect;
	//cs.deCompressYarn(allPlaneIntersect, ellipses, deCompressPlaneIntersect);
	//allPlaneIntersect = deCompressPlaneIntersect;

	//write ellipses to file for testing
	// NOTE: don't rewrite for the procedural case, use simulated compress parameters to validate the compression
	FILE *fout;
	if (fopen_s(&fout, "../data/orientation.txt", "wt") == 0 ) {
		for (int i = 0; i < ellipses.size(); ++i) {
			if (i == 54 || i == 59)
				std::cout << i << "  " << ellipses[i].angle << std::endl;

			fprintf_s(fout, "%.4f %.4f \n", ellipses[i].center.x, ellipses[i].center.y);
			fprintf_s(fout, "%.4f %.4f %.4f \n", ellipses[i].longR, ellipses[i].shortR, ellipses[i].angle);
			//fprintf_s(fout, "%.4f %.4f %.4f \n", 0.0286676, 0.0286676, ellipses[i].angle);
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}

	//write the 2D intersections for tetsting
	
	if (fopen_s(&fout, "../data/allCrossSection2D.txt", "wt") == 0) {
		fprintf_s(fout, "plane_num: %d \n", allPlaneIntersect.size());
		fprintf_s(fout, "ply_num: %d \n", allPlaneIntersect[0].size());
		fprintf_s(fout, "\n");
		for (int i = 0; i < allPlaneIntersect.size(); ++i) { //number of planes
			for (int p = 0; p < allPlaneIntersect[i].size(); ++p) { //number of plys
				fprintf_s(fout, "ply_fiber_num: %d \n", allPlaneIntersect[i][p].size());
				vec2f plyCenter(0.f);
				for (int j = 0; j < allPlaneIntersect[i][p].size(); ++j) { //number of intersections
					plyCenter += allPlaneIntersect[i][p][j];
				}
				plyCenter /= allPlaneIntersect[i][p].size();
				//fprintf_s(fout, "plyCenter: %.4lf %.4lf \n", plyCenter.x, plyCenter.y); // ******* simulated **********
				fprintf_s(fout, "plyCenter: %.4lf %.4lf \n", allPlaneIntersect[i][p][0].x, allPlaneIntersect[i][p][0].y); // ****** Procedural ********
				for (int j = 0; j < allPlaneIntersect[i][p].size(); ++j) { //number of intersections
					fprintf_s(fout, "%.4f %.4f \n", allPlaneIntersect[i][p][j].x, allPlaneIntersect[i][p][j].y);
				}
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
}

void ply_centers_test() {
	//const char* yarnfile = "genYarn.txt"; //For procedural yarn
	const char* yarnfile = "frame00029_scaled.txt"; //For simulated yarn
	const char* curvefile = "frame00029_avg.txt";
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);

	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);

	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);
	
	// Decompress simulated yarn
	//std::vector<yarnIntersect2D> deCompressPlaneIntersect;
	//cs.deCompressYarn(allPlaneIntersect, ellipses, deCompressPlaneIntersect);
	//allPlaneIntersect = deCompressPlaneIntersect;

	//transform from e1-e2 space to x-y space
	std::vector<yarnIntersect2D> allPlaneIntersect_world;
	cs.transferLocal2XY(allPlaneIntersect, allPlaneIntersect_world);

	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect_world, ellipses);
	allPlaneIntersect = allPlaneIntersect_world;

	//extract ply-centers helix parameter
	//cs.extractPlyTwist(allPlaneIntersect, "plyCenter.txt"); //todo

	FILE *fout1;
	//write ellipses to file for testing
	if (fopen_s(&fout1, "../data/orientation.txt", "wt") == 0) {
		const int ignorPlanes = 0.1 * allPlaneIntersect.size(); // crop the first and last 10% of the yarn
		for (int i = ignorPlanes; i < ellipses.size()- ignorPlanes; ++i) {
			fprintf_s(fout1, "%.4f %.4f \n", ellipses[i].center.x, ellipses[i].center.y);
			fprintf_s(fout1, "%.4f %.4f %.4f \n", ellipses[i].longR, ellipses[i].shortR, ellipses[i].angle);
			fprintf_s(fout1, "\n");
		}
		fclose(fout1);
	}

	FILE *fout;
	// write the plycenters
	if (fopen_s(&fout, "../data/allCrossSection2D.txt", "wt") == 0) {
		const int ignorPlanes = 0.1 * allPlaneIntersect.size(); // crop the first and last 10% of the yarn
		/* in case of visualizing ellipses and points together, make sure they are in same amount!!! */
		fprintf_s(fout, "plane_num: %d \n", allPlaneIntersect.size() - 2*ignorPlanes);
		fprintf_s(fout, "ply_num: %d \n", allPlaneIntersect[0].size());
		fprintf_s(fout, "\n");
		
		for (int i = ignorPlanes; i < allPlaneIntersect.size() - ignorPlanes; ++i) { //number of planes
			for (int p = 0; p < allPlaneIntersect[i].size(); ++p) { //number of plys
				fprintf_s(fout, "ply_fiber_num: %d \n", allPlaneIntersect[i][p].size());
				vec2f plyCenter(0.f);
				for (int j = 0; j < allPlaneIntersect[i][p].size(); ++j) { //number of intersections
					plyCenter += allPlaneIntersect[i][p][j];
				}
				plyCenter /= allPlaneIntersect[i][p].size();
				fprintf_s(fout, "plyCenter: %.4lf %.4lf \n", plyCenter.x, plyCenter.y);
				
				for (int j = 0; j < allPlaneIntersect[i][p].size(); ++j) { //number of intersections
					fprintf_s(fout, "%.4f %.4f \n", allPlaneIntersect[i][p][j].x, allPlaneIntersect[i][p][j].y);
				}
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
}					


void extractNormals()
{
	const char* yarnfile = "frame00029_scaled.txt"; //For simulated yarn
	const char* curvefile = "frame00029_avg.txt";
	CrossSection cs(yarnfile, curvefile, 2, 1526, 100);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);

	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect, ellipses);


	//extract spline normals
	std::vector<vec3f> normals;
	cs.extractNormals(ellipses, normals, "../data/junk_norm.txt");
	FILE *fout;
	if (fopen_s(&fout, "../data/normals.txt", "wt") == 0) {
		for (int i = 0; i < ellipses.size(); ++i) {
			fprintf_s(fout, "%.4f %.4f %.4f \n", normals[i].x, normals[i].y, normals[i].z);
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
	// Testing entire curve
	{
		HermiteCurve curve = cs.get_curve();
		const int n = 1526; //normals.size()
		Eigen::Vector3d vtx[n], tang[n], norm[n];
		//curve.output(n, vtx, tang, norm);
		FILE *fout;
		if (fopen_s(&fout, "../data/junk_multiple.txt", "wt") == 0) {
			for (int i = 0; i < n; ++i) {
				Plane plane;
				cs.get_plane(i, plane);
				if (i % 10 == 0) {
					fprintf_s(fout, "%.4lf %.4lf %.4lf ", plane.point.x, plane.point.y, plane.point.z); //curve point
					fprintf_s(fout, "%.4lf %.4lf %.4lf ", plane.n.x, plane.n.y, plane.n.z); //curve tang
					fprintf_s(fout, "%.4lf %.4lf %.4lf\n", normals[i].x, normals[i].y, normals[i].z); //curve normal: read from file
				}
			}
			fclose(fout);
		}
	}
}
