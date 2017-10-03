#include "Fiber.h"

int main() {
	
	Fiber::Yarn yarn;

	yarn.parse("config.txt");
	yarn.yarn_simulate();
	//yarn.compress_yarn("compress.txt");
	yarn.curve_yarn("splines.txt");
	yarn.write_yarn("gen_yarn.txt");

	return 0;
}