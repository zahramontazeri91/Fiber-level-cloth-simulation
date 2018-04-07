#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"
#include <string>


int main(int argc, const char **argv) {

	const char* yarnfile1 = "genYarn_ref.txt";
	//const char* configfile = "configFiles/yarn4/config_step2.txt";
	//const char* configfile = "configFiles/yarn11/config_step2.txt";
	const char* configfile = "configFiles/yarn9/config_step2.txt";
	std::ifstream fin0(configfile);
	assert(fin0.is_open() && "config file wasn't found!\n");
	Fiber::Yarn yarn;
	yarn.parse(configfile);

	yarn.setStepNum(300);
	//yarn.setStepNum(250);
	
	yarn.yarn_simulate();
	yarn.write_yarn(yarnfile1);
	//return 0;

	int phase = 1;

	switch (phase) {
		case 1: {		

			const int window_size = 9;
			/**************** RUN ALL ****************/
			int yarn0 = 0;
			int yarn1 = 1;
			int skipFactor = 500;
			int frame0 = 8000 / skipFactor ;
			int isTrain = 1;
			float trimPercent = 0.1;

			std::cout << " ********************** spacing0.5x *********************** \n" << std::endl;
			std::string dataset = "pattern/yarn11/spacing0.5x/10";
			int frame1 = 14000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing0.5x/00011";
			//frame1 = 16000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing0.5x/10100";
			//frame1 = 15000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing0.5x/11110";
			//frame1 = 15000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//
			//std::cout << " ********************** spacing1.0x *********************** \n" << std::endl;
			//dataset = "pattern/yarn11/spacing1.0x/10";
			//frame1 = 14500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			dataset = "pattern/yarn9/spacing1.0x/00011"; /****************///////*****************/
			frame1 = 17000 / skipFactor + 1;
			full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing1.0x/10100";
			//frame1 = 15500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing1.0x/11110";
			//frame1 = 16000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//
			//std::cout << " ********************** spacing1.5x *********************** \n" << std::endl;
			//dataset = "pattern/yarn11/spacing1.5x/10";
			//frame1 = 15000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing1.5x/00011";
			//frame1 = 17500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing1.5x/10100";
			//frame1 = 16000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing1.5x/11110";
			//frame1 = 16500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);

			////std::cout << " ********************** spacing2.0x *********************** \n" << std::endl;
			//dataset = "pattern/yarn11/spacing2.0x/10";
			//frame1 = 16000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing2.0x/00011";
			//frame1 = 18000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing2.0x/10100";
			//frame1 = 17000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing2.0x/11110";
			//frame1 = 17500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);

			//std::cout << " ********************** spacing2.5x *********************** \n" << std::endl;
			//dataset = "pattern/yarn11/spacing2.5x/10";
			//frame1 = 16500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing2.5x/00011";
			//frame1 = 18500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing2.5x/10100";
			//frame1 = 17500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing2.5x/11110";
			//frame1 = 18500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);

			//std::cout << " ********************** spacing3.0x *********************** \n" << std::endl;
			//dataset = "pattern/yarn11/spacing3.0x/10";
			//frame1 = 16500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing3.0x/00011";
			//frame1 = 19000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing3.0x/10100";
			//frame1 = 20000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			//dataset = "pattern/yarn11/spacing3.0x/11110";
			//frame1 = 19000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);

			std::cout << " ********************** woven *********************** \n" << std::endl;
			yarn0 = 10;
			yarn1 = 11;
			skipFactor = 500;
			isTrain = 0;
			dataset = "woven/yarn4/spacing1.0x/00011";
			frame0 = 0 / skipFactor;
			frame1 = 8000 / skipFactor + 1;
			trimPercent = 0.0;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);

			std::cout << " ********************** twist *********************** \n" << std::endl;
			yarn0 = 10;
			yarn1 = 11;
			skipFactor = 2000;
			isTrain = 0;
			dataset = "twist/yarn4/damp2_500";
			frame0 = 0 / skipFactor;
			frame1 = 99500 / skipFactor + 1;
			trimPercent = 0.0;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);

			std::cout << " ********************** stretch *********************** \n" << std::endl;
			yarn0 = 0;
			yarn1 = 46;
			skipFactor = 500;
			isTrain = 0;
			dataset = "stretch/yarn4/stretch";
			frame0 = 10000 / skipFactor;
			frame1 = 74500 / skipFactor + 1;
			trimPercent = 0.0;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset, isTrain, window_size, trimPercent);
			break;
		}
		case 2: {
			/**************** RUN ALL ****************/
			int yarn0 = 0;
			int yarn1 = 1;
			int skipFactor = 500;		

			std::cout << " ********************** spacing0.5x *********************** \n" << std::endl;
			std::string dataset = "pattern/yarn11/spacing0.5x/10";
			int frame0 = 12000 / skipFactor;
			int frame1 = 14000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing0.5x/00011";
			//frame0 = 14000 / skipFactor;
			//frame1 = 16000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing0.5x/10100";
			//frame0 = 13000 / skipFactor;
			//frame1 = 15000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing0.5x/11110";
			//frame0 = 13000 / skipFactor;
			//frame1 = 15000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//
			//std::cout << " ********************** spacing1.0x *********************** \n" << std::endl;
			//dataset = "pattern/yarn11/spacing1.0x/10";
			//frame0 = 12500 / skipFactor;
			//frame1 = 14500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing1.0x/00011";
			//frame0 = 15000 / skipFactor; 
			//frame1 = 17000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing1.0x/10100";
			//frame0 = 13500 / skipFactor;
			//frame1 = 15500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing1.0x/11110";
			//frame0 = 14000 / skipFactor;
			//frame1 = 16000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);

			//std::cout << " ********************** spacing1.5x *********************** \n"<< std::endl;
			//dataset = "pattern/yarn11/spacing1.5x/10";
			//frame0 = 13000 / skipFactor;
			//frame1 = 15000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing1.5x/00011";
			//frame0 = 15500 / skipFactor;
			//frame1 = 17500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing1.5x/10100";
			//frame0 = 14000 / skipFactor;
			//frame1 = 16000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing1.5x/11110";
			//frame0 = 14500 / skipFactor;
			//frame1 = 16500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);

			//std::cout << " ********************** spacing2.0x *********************** \n" << std::endl;
			//dataset = "pattern/yarn11/spacing2.0x/10";
			//frame0 = 14000 / skipFactor;
			//frame1 = 16000 / skipFactor + 1;
			////step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing2.0x/00011";
			//frame0 = 16000 / skipFactor;
			//frame1 = 18000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing2.0x/10100";
			//frame0 = 15000 / skipFactor;
			//frame1 = 17000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing2.0x/11110";
			//frame0 = 15500 / skipFactor;
			//frame1 = 17500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);

			//std::cout << " ********************** spacing2.5x *********************** \n" << std::endl;
			//dataset = "pattern/yarn11/spacing2.5x/10";
			//frame0 = 14500 / skipFactor; 
			//frame1 = 16500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing2.5x/00011";
			//frame0 = 16500 / skipFactor;
			//frame1 = 18500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing2.5x/10100";
			//frame0 = 15500 / skipFactor;
			//frame1 = 17500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing2.5x/11110";
			//frame0 = 16500 / skipFactor;
			//frame1 = 18500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);

			//std::cout << " ********************** spacing3.0x *********************** \n" << std::endl;
			//dataset = "pattern/yarn11/spacing3.0x/10";
			//frame0 = 14500 / skipFactor;
			//frame1 = 16500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing3.0x/00011";
			//frame0 = 17000 / skipFactor;
			//frame1 = 19000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing3.0x/10100";
			//frame0 = 18000 / skipFactor;
			//frame1 = 20000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);
			//dataset = "pattern/yarn11/spacing3.0x/11110";
			//frame0 = 17000 / skipFactor;
			//frame1 = 19000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);

			std::cout << " ********************** woven *********************** \n" << std::endl;
			yarn0 = 10;
			yarn1 = 11;
			skipFactor = 500;
			dataset = "woven/yarn4/spacing1.0x/00011";
			frame0 = 0 / skipFactor;
			frame1 = 8000 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);


			std::cout << " ********************** twist *********************** \n" << std::endl;
			yarn0 = 10;
			yarn1 = 11;;
			skipFactor = 2000;
			dataset = "twist/yarn4/damp2_500";
			frame0 = 0 / skipFactor;
			frame1 = 99500 / skipFactor + 1;
			//step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);

			std::cout << " ********************** stretch *********************** \n" << std::endl;
			yarn0 = 0;
			yarn1 = 46;
			skipFactor = 500;
			dataset = "stretch/yarn4/stretch";
			frame0 = 34000 / skipFactor;
			frame1 = 74500 / skipFactor + 1;
			step4_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarn0, yarn1, dataset);

			break;
		}
	}

	//	std::system("pause"); //add breakpoint instead 	return 0;
}