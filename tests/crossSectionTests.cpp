#include "hermiteTests.h"
#include "crossSectionTests.h"
#include "../crossSection.h"

void linePlaneIntersection_test() {
	const char* yarnfile = "D:/sandbox/fiberSimulation/yarn_generation_project/results/output00001.txt";
	const char* curvefile = "avg00029.txt"; // TO DO: why avg is wrong TODO: cleanup input files
	CrossSection cs(yarnfile, curvefile, 10,5);

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
	const char* yarnfile = "D:/sandbox/fiberSimulation/yarn_generation_project/results/output00001.txt";
	const char* curvefile = "avg00029.txt"; 
	CrossSection cs(yarnfile, curvefile, 10,5);

	Plane plane;
	plane.normal = vec3f(0,0,1);
	plane.point = vec3f(0, 0, 0);
	std::vector<vec3f> itsList;
	bool b = cs.yarnPlaneIntersection(plane, itsList);

	FILE *fout;
	if (fopen_s(&fout, "crossSection.txt", "wt") == 0) {
		fprintf_s(fout, "%d \n", itsList.size());
		// First write the yarn-center
		fprintf_s(fout, "%.4lf %.4lf %.4lf \n", plane.point[0], plane.point[1], plane.point[0]);
		for (int i = 0; i < itsList.size(); ++i) {
			fprintf_s(fout, "%.4lf %.4lf %.4lf \n", itsList[i][0], itsList[i][1], itsList[i][2]);
		}
		fclose(fout);
	}
}

void bildPlanes_test() {
	const char* yarnfile = "D:/sandbox/fiberSimulation/yarn_generation_project/results/output00001.txt";
	const char* curvefile = "avg00029.txt";
	CrossSection cs(yarnfile, curvefile, 10, 5);

	Plane plane;
	plane.normal = vec3f(0, 0, 1);
	plane.point = vec3f(0, 0, 0);
	std::vector<vec3f> itsList;
	bool b = cs.yarnPlaneIntersection(plane, itsList);
}