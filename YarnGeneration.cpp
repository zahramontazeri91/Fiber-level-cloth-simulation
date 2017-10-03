#include "Fiber.h"
#include "tests/hermiteTests.h"

int main() {
#if 0	
	Fiber::Yarn yarn;

	yarn.parse("config.txt");
	yarn.yarn_simulate();
	yarn.compress_yarn("compress.txt");
	yarn.curve_yarn("splines.txt");
	yarn.write_yarn("gen_yarn.txt");
#else
    hermiteTest1();
    hermiteTest2();
#endif

	return 0;
}
