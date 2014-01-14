#include "Sequences.h"
#include "FeaturesProfile.h"
#include "Profile.h"
#include "txtProc.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
namespace po = boost::program_options;
int main(int argc, char *argv[]){
	int gapPen,codonLength, phosphScore,domainScore;
	double identityCutOff, gapExt;
	bool weightsModeOn;
	std::string filename,verboseMode,outputPrefix;
	po::options_description desc("Allowed options");
	desc.add_options()
    		("help,h", "produce help message")
    		("input,i", po::value<std::string>(&filename), "input file name")
    		("output,o", po::value<std::string>(&outputPrefix), "output file prefix")
		("gap_penalty,g", po::value<int>(&gapPen), "gap opening penalty")
		("gap_extension,e",po::value<double>(&gapExt)->default_value(1.),"gap extension penalty")
		("codon_length,c", po::value<int>(&codonLength)->implicit_value(4)->default_value(1),"codon length")
		("phosph,p", po::value<int>(&phosphScore)->default_value(15),"score for aligning phosphoryated residues")
		("domain,d", po::value<int>(&domainScore)->default_value(2),"score for aligning domains")
		("weights,w", po::value<bool>(&weightsModeOn)->implicit_value(true)->default_value(false),"all sequences contribute to the profile with weights(=similarity)")
		("identity", po::value<double>(&identityCutOff)->default_value(0.8),"identity cut off for sequences included in profile")
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
			//Profile prf;													//this prf will be useful for next rounds of alignments
			//FeaturesProfile fprf;
			std::cout << weightsModeOn << std::endl;
			std::vector<std::vector<double> > prf, fprf;
			std::vector<std::string> multipleAlignment(rawSequences.performMSAencoded(&prf,&fprf,gapPen,gapExt,verboseMode,weightsModeOn));	//create multiple sequence alignment	

			std::vector<std::vector<std::vector<std::string> > > encSeq = rawSequences.getEncodedSequences();
			std::string outputPrefix1 = "crap1_step1";
			//txtProc::writeAlignmentWithoutCodeToFile(multipleAlignment,encSeq,outputPrefix1);						//write multiple alignment to a fileA
			txtProc::writeAlignmentToFile(multipleAlignment,encSeq,outputPrefix1);						//write multiple alignment to a fileA
			Profile prof(prf);
			FeaturesProfile fprof(fprf);
			std::vector<std::string> multipleAlignment2ndRound;
			for (int i = 10; i > 0; i--){
				double cutoff = 1/i;
				multipleAlignment2ndRound=rawSequences.performMSAnextRound(&prof,&fprof,gapPen,gapExt,verboseMode,weightsModeOn,cutoff);
			}
			//txtProc::writeAlignmentWithoutCodeToFile(multipleAlignment2ndRound,encSeq,outputPrefix);						//write multiple alignment to a fileA
			
			multipleAlignment2ndRound=rawSequences.performMSAnextRound(&prof,&fprof,gapPen,gapExt,verboseMode,weightsModeOn,0);
			txtProc::writeAlignmentToFile(multipleAlignment2ndRound,encSeq,"profile0cutoff");						//write multiple alignment to a fileA
			prof.printProfile(0,1);
			txtProc::writeAlignmentToFile(multipleAlignment2ndRound,encSeq,outputPrefix);						//write multiple alignment to a fileA
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
