#include <fstream>
#include "Fiber.h"
#include "crossSection.h"
#include "fitting.h"

#include "tests/hermiteTests.h"
#include "tests/crossSectionTests.h"
#include <string>


int main(int argc, const char **argv) {

	const char* yarnfile1 = "genYarn_ref.txt";
	const char* configfile = "config_300.txt";
	std::ifstream fin0(configfile);
	assert(fin0.is_open() && "config file wasn't found!\n");
	Fiber::Yarn yarn;
	yarn.parse(configfile);

	yarn.setStepNum(300);
	
	yarn.yarn_simulate();
	yarn.write_yarn(yarnfile1);

	int phase = 1;

	switch (phase) {
		case 1: {
			const int yarnNum = 1;
			const int skipFactor = 500;

			/**************** RUN ALL ****************/
			const int frame0 = 17000 / skipFactor ;
			const int isTrain = 1;
			const int window_size = 50;
			const float trimPercent = 0.1;

			std::cout << " ********************** spacing0.5x *********************** \n" << std::endl;
			std::string dataset = "yarn4/spacing0.5x/10";
			int frame1 = 14000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			dataset = "yarn4/spacing0.5x/00011";
			frame1 = 16000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			dataset = "yarn4/spacing0.5x/10100";
			frame1 = 15000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			dataset = "yarn4/spacing0.5x/11110";
			frame1 = 15000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			
			std::cout << " ********************** spacing1.0x *********************** \n" << std::endl;
			dataset = "yarn4/spacing1.0x/10";
			frame1 = 14500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			dataset = "yarn4/spacing1.0x/00011";
			frame1 = 17000 / skipFactor + 1;
			full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			dataset = "yarn4/spacing1.0x/10100";
			frame1 = 15500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			dataset = "yarn4/spacing1.0x/11110";
			frame1 = 16000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			
			std::cout << " ********************** spacing1.5x *********************** \n" << std::endl;
			dataset = "yarn4/spacing1.5x/10";
			frame1 = 15000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			dataset = "yarn4/spacing1.5x/00011";
			frame1 = 17500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			dataset = "yarn4/spacing1.5x/10100";
			frame1 = 16000 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);
			dataset = "yarn4/spacing1.5x/11110";
			frame1 = 16500 / skipFactor + 1;
			//full_pipeline(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset, isTrain, window_size, trimPercent);

			break;
		}
		case 2: {
			const int yarnNum = 1;
			const int skipFactor = 500;
			/**************** RUN ALL ****************/
			const int frame0 = 17000 / skipFactor;		

			std::cout << " ********************** spacing0.5x *********************** \n" << std::endl;
			std::string dataset = "yarn4/spacing0.5x/10";
			int frame1 = 14000 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "yarn4/spacing0.5x/00011";
			frame1 = 16000 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "yarn4/spacing0.5x/10100";
			frame1 = 15000 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "yarn4/spacing0.5x/11110";
			frame1 = 15000 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			
			std::cout << " ********************** spacing1.0x *********************** \n" << std::endl;
			dataset = "yarn4/spacing1.0x/10";
			frame1 = 14500 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "yarn4/spacing1.0x/00011";
			frame1 = 17000 / skipFactor + 1;
			step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "yarn4/spacing1.0x/10100";
			frame1 = 15500 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "yarn4/spacing1.0x/11110";
			frame1 = 16000 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);

			std::cout << " ********************** spacing1.5x *********************** \n"<< std::endl;
			dataset = "yarn4/spacing1.5x/10";
			frame1 = 15000 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "yarn4/spacing1.5x/00011";
			frame1 = 17500 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "yarn4/spacing1.5x/10100";
			frame1 = 16000 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);
			dataset = "yarn4/spacing1.5x/11110";
			frame1 = 16500 / skipFactor + 1;
			//step3_NN_output(yarnfile1, configfile, yarn, skipFactor, frame0, frame1, yarnNum, dataset);

			break;
		}
	}

	//	std::system("pause"); //add breakpoint instead

	return 0;
}