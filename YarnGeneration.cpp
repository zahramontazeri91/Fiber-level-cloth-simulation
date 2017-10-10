#include <fstream>
#include "Fiber.h"
#include "tests/hermiteTests.h"
#include "crossSection.h"

int main(int argc, const char **argv) {
#if 0
    
	//Fiber::Yarn yarn;
	//yarn.parse("config_v2.txt");
	//yarn.yarn_simulate();
	//yarn.compress_yarn("compress.txt");
	//yarn.curve_yarn("curves.txt");
	//yarn.write_yarn("gen_yarn.txt");
    

    if ( argc != 3 ) {
        printf("Usage: YarnGeneration [task file] [output file]\n");
		return 1;
    }

    std::ifstream fin(argv[1]);
    if ( fin.is_open() ) {
        std::cout << "Using task file: \"" << argv[1] << "\"." << std::endl;

        Fiber::Yarn yarn;

        std::string command, fname;
        fin >> fname;
        yarn.parse(fname.c_str());
        yarn.yarn_simulate();

        while ( fin >> command >> fname )
            if ( command == "COMPRESS" )
                yarn.compress_yarn(fname.c_str());
            else if ( command == "CURVE" )
                yarn.curve_yarn(fname.c_str());

        yarn.write_yarn(argv[2]);
    }


#else
    //hermiteTest1();
    //hermiteTest2();

	const char* yarnfile = "D:/sandbox/fiberSimulation/yarn_generation_project/results/output00001.txt";
	const char* curvefile = "curves.txt"; // TODO: D: / sandbox / fiberSimulation / hairs / avg_hair / frame00029_yarn.txt";
	CrossSection cs (yarnfile, curvefile, 10);
	while (1);

#endif

	return 0;
}
