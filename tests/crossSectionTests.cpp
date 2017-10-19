#include "hermiteTests.h"
#include "crossSectionTests.h"
#include "../crossSection.h"

using namespace std;

void linePlaneIntersection_test() {
	const char* yarnfile = "../data/gen_yarn_f1.txt";
	const char* curvefile = "../data/frame00029_avg.txt";
	CrossSection cs(yarnfile, curvefile, 2, 10000, 100);

	Plane plane;
	plane.normal = vec3f(1, 0, 0);
	plane.point = vec3f(0, 0, 0);
	vec3f its(0.f);
	vec3f start(1, 0, 0);
	vec3f end(-1, 0, 0);
	bool a = cs.linePlaneIntersection(start, end, plane, its);
	std::cout << " Intersection: " << a << "  " << its << std::endl;
}

void yarnPlaneIntersection_test() {
	const char* yarnfile = "../data/gen_yarn_f1.txt";
	const char* curvefile = "../data/frame00029_avg.txt"; 
	CrossSection cs(yarnfile, curvefile, 2, 10000, 100);

	Plane plane;
	cs.get_plane(1, plane);

	//Plane plane;
	//plane.normal = vec3f(-0.0316426, -0.304221,  0.952076);
	//plane.point = vec3f(0,0, plane0.point.z);

	std::cout << plane.normal.x << "  " << plane.normal.y << "  " << plane.normal.z << std::endl;
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
	const char* yarnfile = "../data/gen_yarn_f1.txt";
	const char* curvefile = "../data/frame00029_avg.txt";
	const int num_planes = 1000;
	CrossSection cs(yarnfile, curvefile, 2, 10000, 100);
	FILE *fout;
	if (fopen_s(&fout, "../data/test_planes.txt", "wt") == 0) {
		fprintf_s(fout, "%d \n", num_planes);
		for (int i = 0; i < num_planes; i++) {
			Plane plane;
			cs.get_plane(i, plane);
			fprintf_s(fout, "%.4lf %.4lf %.4lf \n", plane.normal.x, plane.normal.y, plane.normal.z);
		}
		fclose(fout);
	}
}

void allPlanesIntersections_test() {
	const char* yarnfile = "../data/gen_yarn_f1.txt"; //For procedural yarn
	//const char* yarnfile = "../data/frame00001_scaled.txt"; //For simulated yarn
	const char* curvefile = "../data/frame00001_avg.txt"; 
	CrossSection cs(yarnfile, curvefile, 2, 10000, 100);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	//std::cout << "intersections lists size: " << itsLists.size() << std::endl;

	cs.write_PlanesIntersections3D("../data/allCrossSection.txt",itsLists);
}

void project2Plane_test() {
	const char* yarnfile = "../data/frame00001_scaled.txt"; //For simulated yarn
	const char* curvefile = "../data/frame00029_avg.txt"; 
	CrossSection cs(yarnfile, curvefile, 2, 10000, 100);
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
	const char* yarnfile = "../data/gen_yarn_f1.txt"; //For procedural yarn
	//const char* yarnfile = "../data/frame00001_scaled.txt"; //For simulated yarn
	const char* curvefile = "../data/frame00001_avg.txt"; 
	CrossSection cs(yarnfile, curvefile, 2, 10000, 100);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	std::cout << "intersections lists size: " << itsLists.size() << std::endl;

	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);
	std::cout << "number of intersected planes: " << allPlaneIntersect.size() << std::endl;

	//write the 2D intersections for tetsting
	FILE *fout;
	if (fopen_s(&fout, "../data/allCrossSection2D.txt", "wt") == 0) {
		fprintf_s(fout, "plane_num: %d \n", allPlaneIntersect.size());
		fprintf_s(fout, "ply_num: %d \n", allPlaneIntersect[0].size());
		fprintf_s(fout, "\n");
		for (int i = 0; i < allPlaneIntersect.size(); ++i) { //number of planes
			//Plane plane;
			//cs.get_plane(i, plane);
			//fprintf_s(fout, "center: %.4lf %.4lf %.4lf \n", plane.point[0], plane.point[1], plane.point[2]); //center is always 0,0
			for (int p = 0; p < allPlaneIntersect[i].size(); ++p) { //number of plys
				fprintf_s(fout, "ply_fiber_num: %d \n", allPlaneIntersect[i][p].size());
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
	const char* yarnfile = "../data/gen_yarn_f1.txt"; //For procedural yarn
	//const char* yarnfile = "../data/frame00001_scaled.txt"; //For simulated yarn
	const char* curvefile = "../data/frame00001_avg.txt";
	CrossSection cs(yarnfile, curvefile, 2, 10000, 100);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);

	Ellipse e;
	cs.fitEllipse(allPlaneIntersect[5], e);
	FILE *fout;
	if (fopen_s(&fout, "../data/pca_test.txt", "wt") == 0) {
		for (int p = 0; p < allPlaneIntersect[5].size(); ++p) {
			fprintf_s(fout, "%d \n", allPlaneIntersect[5][p].size());
			for (int i = 0; i < allPlaneIntersect[5][p].size(); ++i) {
				fprintf_s(fout, "%.6f %.6f \n", allPlaneIntersect[5][p][i].x, allPlaneIntersect[5][p][i].y);
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
	//const char* yarnfile = "../data/gen_yarn_f1.txt"; //For procedural yarn
	const char* yarnfile = "../data/frame00029_scaled.txt"; //For simulated yarn
	const char* curvefile = "../data/frame00029_avg.txt";
	CrossSection cs(yarnfile, curvefile, 2, 10000, 100);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	std::vector<yarnIntersect2D> allPlaneIntersect;
	cs.PlanesIntersections2D(itsLists, allPlaneIntersect);

	std::vector<Ellipse> ellipses;
	cs.extractCompressParam(allPlaneIntersect, ellipses, "compress.txt");

	//write ellipses to file for testing
	FILE *fout;
	if (fopen_s(&fout, "../data/orientation.txt", "wt") == 0 ) {
		for (int i = 0; i < ellipses.size(); ++i) {
			fprintf_s(fout, "%.4f %.4f \n", ellipses[i].center.x, ellipses[i].center.y);
			fprintf_s(fout, "%.4f %.4f %.4f \n", ellipses[i].longR, ellipses[i].shortR, ellipses[i].angle);
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}


	//write the 2D intersections for tetsting
	//FILE *fout;
	if (fopen_s(&fout, "../data/allCrossSection2D.txt", "wt") == 0) {
		fprintf_s(fout, "plane_num: %d \n", allPlaneIntersect.size());
		fprintf_s(fout, "ply_num: %d \n", allPlaneIntersect[0].size());
		fprintf_s(fout, "\n");
		for (int i = 0; i < allPlaneIntersect.size(); ++i) { //number of planes
			for (int p = 0; p < allPlaneIntersect[i].size(); ++p) { //number of plys
				fprintf_s(fout, "ply_fiber_num: %d \n", allPlaneIntersect[i][p].size());
				for (int j = 0; j < allPlaneIntersect[i][p].size(); ++j) { //number of intersections
					fprintf_s(fout, "%.4f %.4f \n", allPlaneIntersect[i][p][j].x, allPlaneIntersect[i][p][j].y);
				}
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}

}

void compress_yarn_test() {
	Fiber::Yarn yarn;

	yarn.parse("config.txt");
	yarn.yarn_simulate();
	yarn.compress_yarn("compress.txt");
	yarn.curve_yarn("frame00029_avg.txt");
	yarn.write_yarn("genYarn.txt");
}