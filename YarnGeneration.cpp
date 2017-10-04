#include "Fiber.h"
#include "tests/hermiteTests.h"

int main() {
#if 1
	Fiber::Yarn yarn;

	yarn.parse("config.txt");
	yarn.yarn_simulate();
	yarn.compress_yarn("compress.txt");
	yarn.curve_yarn("curves.txt");
	yarn.write_yarn("gen_yarn.txt");
#else
    hermiteTest1();
    hermiteTest2();
#endif

	return 0;
}
