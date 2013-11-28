#include "Sequences.h"
#include "Profile.h"
#include "txtProc.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
namespace po = boost::program_options;
int main(int argc, char *argv[]){
	int gapPen;
	std::string filename;
	std::string verboseMode;
	po::options_description desc("Allowed options");
	desc.add_options()
    		("help", "produce help message")
    		("i", po::value<std::string>(&filename), "input file name")
		("g", po::value<int>(&gapPen), "gap opening penalty")
		("v",po::value<std::string>(&verboseMode)->implicit_value("1")->default_value("0"),"verbose mode")
	;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if (vm.count("help")) {
    		std::cout << desc << "\n";
    		return 1;
	}
	Sequences rawSequences(txtProc::processFASTA(filename));			//read data from file
	Profile prf;										//this prf will be useful for next rounds of alignments
	std::vector<std::string> multipleAlignment(rawSequences.performMSA(&prf,gapPen,verboseMode));	//create multiple sequence alignment	
	txtProc::writeAlignmentToFile(multipleAlignment,rawSequences.getSequences(),filename);	//write multiple alignment to a file
}
