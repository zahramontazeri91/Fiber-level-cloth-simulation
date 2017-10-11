#include "hermiteTests.h"
#include "crossSectionTests.h"
#include "../crossSection.h"

void linePlaneIntersection_test() {
	const char* yarnfile = "D:/sandbox/fiberSimulation/yarn_generation_project/results/output00029.txt";
	const char* curvefile = "avg00029.txt"; // TO DO: why avg is wrong TODO: cleanup input files
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
	const char* yarnfile = "D:/sandbox/fiberSimulation/yarn_generation_project/results/output00029.txt";
	const char* curvefile = "avg00029.txt"; 
	CrossSection cs(yarnfile,2, curvefile, 10,5);

	Plane plane;
	cs.get_plane(1, plane);

	//Plane plane;
	//plane.normal = vec3f(-0.0316426, -0.304221,  0.952076);
	//plane.point = vec3f(0,0, plane0.point.z);

	std::cout << plane.normal.x << "  " << plane.normal.y << "  " << plane.normal.z << std::endl;
	std::cout << plane.point.x << "  " << plane.point.y << "  " << plane.point.z << std::endl;

	std::vector<vec3f> itsList;
	bool b = cs.yarnPlaneIntersection(plane, itsList);

	std::cout << itsList.size() << std::endl;

	FILE *fout;
	if (fopen_s(&fout, "crossSection.txt", "wt") == 0) {
		fprintf_s(fout, "%d \n", itsList.size());
		// First write the yarn-center
		fprintf_s(fout, "%.4lf %.4lf %.4lf \n", plane.point[0], plane.point[1], plane.point[2]);
		for (int i = 0; i < itsList.size(); ++i) {
			fprintf_s(fout, "%.4lf %.4lf %.4lf \n", itsList[i][0], itsList[i][1], itsList[i][2]);
		}
		fclose(fout);
	}
}

void bildPlanes_test() {
	const char* yarnfile = "D:/sandbox/fiberSimulation/yarn_generation_project/results/output00029.txt";
	const char* curvefile = "avg00029.txt";
	const int num_planes = 1000;
	CrossSection cs(yarnfile,2, curvefile, 10, num_planes);
	FILE *fout;
	if (fopen_s(&fout, "test_planes.txt", "wt") == 0) {
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
	const char* yarnfile = "D:/sandbox/fiberSimulation/yarn_generation_project/results/output00029.txt";
	const char* curvefile = "avg00029.txt"; // TO DO: why avg is wrong TODO: cleanup input files
	CrossSection cs(yarnfile,2, curvefile, 10, 5);
	std::vector<std::vector<vec3f>> itsLists;
	cs.allPlanesIntersections(itsLists);
	std::cout << "intersections lists size: " << itsLists.size() << std::endl;

	cs.write_PlanesIntersections("allCrossSection.txt",itsLists);
}