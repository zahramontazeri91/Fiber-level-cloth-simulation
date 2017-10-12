#include "hermiteTests.h"
#include "crossSectionTests.h"
#include "../crossSection.h"

void linePlaneIntersection_test() {
	const char* yarnfile = "../data/output00029.txt";
	const char* curvefile = "../data/frame00029_avg.txt"; // TO DO: why avg is wrong TODO: cleanup input files
	CrossSection cs(yarnfile,2, curvefile, 10,5);

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
	const char* yarnfile = "../data/output00029.txt";
	const char* curvefile = "../data/frame00029_avg.txt"; 
	CrossSection cs(yarnfile,2, curvefile, 10,5);

	Plane plane;
	cs.get_plane(1, plane);

	//Plane plane;
	//plane.normal = vec3f(-0.0316426, -0.304221,  0.952076);
	//plane.point = vec3f(0,0, plane0.point.z);

	std::cout << plane.normal.x << "  " << plane.normal.y << "  " << plane.normal.z << std::endl;
	std::cout << plane.point.x << "  " << plane.point.y << "  " << plane.point.z << std::endl;

	yarnIntersect itsList;
	bool b = cs.yarnPlaneIntersection(plane, itsList);

	std::cout << itsList[0].size() << std::endl;

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

void bildPlanes_test() {
	const char* yarnfile = "../data/output00029.txt";
	const char* curvefile = "../data/frame00029_avg.txt";
	const int num_planes = 1000;
	CrossSection cs(yarnfile,2, curvefile, 10, num_planes);
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
	const char* yarnfile = "../data/output00029.txt"; //For procedural yarn
	//const char* yarnfile = "../data/frame00001_hairs.txt"; //For simulated yarn
	const char* curvefile = "../data/frame00001_avg.txt"; // TO DO: why avg is wrong TODO: cleanup input files
	CrossSection cs(yarnfile,2, curvefile, 100, 1000);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	//std::cout << "intersections lists size: " << itsLists.size() << std::endl;

	cs.write_PlanesIntersections3D("../data/allCrossSection.txt",itsLists);
}

void project2Plane_test() {
	const char* yarnfile = "../data/frame00029_hairs.txt"; //For simulated yarn
	const char* curvefile = "../data/frame00029_avg.txt"; // TO DO: why avg is wrong TODO: cleanup input files
	CrossSection cs(yarnfile, 2, curvefile, 100, 1000);
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
	//const char* yarnfile = "../data/output00001.txt"; //For procedural yarn
	const char* yarnfile = "../data/frame00001_hairs.txt"; //For simulated yarn
	const char* curvefile = "../data/frame00001_avg.txt"; // TODO: cleanup input files
	CrossSection cs(yarnfile, 2, curvefile, 100, 1000);
	std::vector<yarnIntersect> itsLists;
	cs.allPlanesIntersections(itsLists);
	std::cout << "intersections lists size: " << itsLists.size() << std::endl;

	cs.write_PlanesIntersections3D("../data/allCrossSection.txt", itsLists);
	cs.write_PlanesIntersections2D("../data/allCrossSection2D.txt", itsLists);
}