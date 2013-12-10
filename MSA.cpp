#include "Sequences.h"
#include "FeaturesProfile.h"
#include "Profile.h"
#include "txtProc.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
namespace po = boost::program_options;
int main(int argc, char *argv[]){
	int gapPen,codonLength;
	std::string filename,verboseMode,outputPrefix;
	po::options_description desc("Allowed options");
	desc.add_options()
    		("help,h", "produce help message")
    		("input,i", po::value<std::string>(&filename), "input file name")
    		("output,o", po::value<std::string>(&outputPrefix), "output file prefix")
		("gap_penalty,g", po::value<int>(&gapPen), "gap opening penalty")
		("codon_length,c", po::value<int>(&codonLength)->implicit_value(4)->default_value(1),"codon length")
		("verbose,v",po::value<std::string>(&verboseMode)->implicit_value("1")->default_value("0"),"verbose mode")
	;	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if (vm.count("help")) {
    		std::cout << desc << std::endl;
    		return 1;
	}
	if (vm.count("input") && vm.count("gap_penalty")){
		if (codonLength > 1) {
			Sequences rawSequences(txtProc::processFASTA(filename,codonLength));
			Profile prf;													//this prf will be useful for next rounds of alignments
			FeaturesProfile fprf;
			std::vector<std::string> multipleAlignment(rawSequences.performMSAencoded(&prf,&fprf,gapPen,verboseMode));	//create multiple sequence alignment	
			std::vector<std::vector<std::vector<std::string> > > encSeq = rawSequences.getEncodedSequences();
			txtProc::writeAlignmentToFile(multipleAlignment,encSeq,outputPrefix);						//write multiple alignment to a fileA
		}
		else {
			Sequences rawSequences(txtProc::processFASTA(filename));							//read data from file
			Profile prf;													//this prf will be useful for next rounds of alignments
			std::vector<std::string> multipleAlignment(rawSequences.performMSA(&prf,gapPen,verboseMode));			//create multiple sequence alignment	
			txtProc::writeAlignmentToFile(multipleAlignment,rawSequences.getSequences(),outputPrefix);				//write multiple alignment to a file
		}
	}
	else{
		std::cout << "input file and/or gap penalty not specified, try --help for more information" << std::endl;
	}
}
