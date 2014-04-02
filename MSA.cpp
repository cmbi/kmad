#include "Sequences.h"
#include "FeaturesProfile.h"
#include "Profile.h"
#include "txtProc.h"
#include "vecUtil.h"
#include "misc.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <string>
namespace po = boost::program_options;
int main(int argc, char *argv[]){
	int codonLength, phosphScore,domainScore;
	double identityCutOff, gapExt, gapPen;
	bool weightsModeOn;
	std::string filename,verboseMode,outputPrefix;
	po::options_description desc("Allowed options");
	desc.add_options()
	    	("help,h", "produce help message")
	    	("input,i", po::value<std::string>(&filename), "input file name")
	    	("output,o", po::value<std::string>(&outputPrefix), "output file prefix")
		("gap_penalty,g", po::value<double>(&gapPen), "gap opening penalty")
		("gap_extension,e",po::value<double>(&gapExt)->default_value(-1.),"gap extension penalty")
		("codon_length,c", po::value<int>(&codonLength)->implicit_value(4)->default_value(1),"codon length")
		("phosph,p", po::value<int>(&phosphScore)->default_value(15),"score for aligning phosphoryated residues")
		("domain,d", po::value<int>(&domainScore)->default_value(3),"score for aligning domains")
		("weights,w", po::value<bool>(&weightsModeOn)->implicit_value(true)->default_value(false),"all sequences contribute to the profile with weights(=similarity)")
		("identity", po::value<double>(&identityCutOff)->default_value(0.8, "0.8"),"identity cut off for sequences included in profile")
		("verbose,v",po::value<std::string>(&verboseMode)->implicit_value("1")->default_value("0"),"verbose mode");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if (vm.count("help")) {
    		std::cout << desc << std::endl;
    		return 1;
	}
	if (vm.count("input") && vm.count("gap_penalty") && misc::checkParameters(codonLength,phosphScore,domainScore,gapExt,gapPen,weightsModeOn,identityCutOff)){
		if (codonLength > 1) {
			std::vector<std::string> motifs_ids;
			std::vector<double> motifs_probs, identities;
			Sequences rawSequences(txtProc::processFASTA(filename,codonLength, &motifs_ids, &motifs_probs));
			Profile prf;
			FeaturesProfile fprf(domainScore,phosphScore);
			//first round of the alignment - all vs 1st
			std::vector<std::string> multipleAlignment(rawSequences.performMSAencoded(prf,fprf,gapPen,gapExt,verboseMode,weightsModeOn,domainScore,phosphScore,codonLength,identities));
			std::vector<std::vector<std::vector<std::string> > > encSeq = rawSequences.getEncodedSequences();
			Profile prof = prf;
			FeaturesProfile fprof = fprf;
			std::vector<std::string> multipleAlignment2ndRound;
			double gapPenDecreasing;
			for (int i = 10; i >= 0; i--){
				double cutoff = i/10;
				//gapPenDecreasing = -3.-4.*i/10;
				//std::cout << gapPenDecreasing << std::endl;
				multipleAlignment2ndRound=rawSequences.performMSAnextRound(prf,fprf,gapPen,gapExt,verboseMode,weightsModeOn,cutoff,domainScore,phosphScore,codonLength,identities);
				//multipleAlignment2ndRound=rawSequences.performMSAnextRound(prf,fprf,gapPenDecreasing,gapExt,verboseMode,weightsModeOn,cutoff,domainScore,phosphScore);
			}
			//txtProc::writeAlignmentWithoutCodeToFile(multipleAlignment2ndRound,encSeq,outputPrefix);						//write multiple alignment to a fileA
			//gapPenDecreasing = -2;
			txtProc::writeVector(prof.getMatrix(),"dist2");
			multipleAlignment2ndRound=rawSequences.performMSAnextRound(prf,fprf,gapPen,gapExt,verboseMode,weightsModeOn,0,domainScore,phosphScore,codonLength, identities);
			//multipleAlignment2ndRound=rawSequences.performMSAnextRound(prf,fprf,gapPenDecreasing,gapExt,verboseMode,weightsModeOn,0,domainScore,phosphScore);
			txtProc::writeAlignmentToFile(multipleAlignment2ndRound,encSeq,outputPrefix);						//write multiple alignment to a fileA
		}
		else {
			Sequences rawSequences(txtProc::processFASTA(filename));							//read data from file
			Profile prf;													//this prf will be useful for next rounds of alignments
			std::vector<std::string> multipleAlignment(rawSequences.performMSA(prf,gapPen,verboseMode));			//create multiple sequence alignment
			txtProc::writeAlignmentToFile(multipleAlignment,rawSequences.getSequences(),outputPrefix);				//write multiple alignment to a file
		}
	}
	else{
		std::cout << "input file and/or gap penalty not specified, try --help for more information" << std::endl;
	}
}
