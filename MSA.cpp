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
	int codonLength, phosphScore,domainScore, motifScore;
	double identityCutOff, gapExt, gapPen, lcr_mod;
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
		("phosph,p", po::value<int>(&phosphScore)->default_value(0),"score for aligning phosphoryated residues")
		("domain,d", po::value<int>(&domainScore)->default_value(0),"score for aligning domains")
		("lcr,l", po::value<double>(&lcr_mod)->default_value(1), "gap penalty modifier inside a low complexity region")
		("motif,m", po::value<int>(&motifScore),"probability multiplier for motifs")
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
	if (vm.count("input") && vm.count("gap_penalty") && misc::checkParameters(codonLength,phosphScore,domainScore,motifScore,gapExt,gapPen,weightsModeOn,identityCutOff)){
		std::vector<std::string> motifs_ids;
		std::vector<double> motifs_probs, identities;
		Sequences rawSequences(txtProc::processFASTA(filename,codonLength, &motifs_ids, &motifs_probs));
		Profile prf;
		FeaturesProfile fprf(domainScore,phosphScore,motifScore,lcr_mod,motifs_ids,motifs_probs);
		//first round of the alignment - all vs 1st
		std::vector<std::string> multipleAlignment(rawSequences.performMSAencoded(prf,fprf,gapPen,gapExt,verboseMode,weightsModeOn,codonLength,identities));
		std::vector<std::vector<std::vector<std::string> > > encSeq = rawSequences.getEncodedSequences();
		std::vector<std::string> multipleAlignment2ndRound;
		double gapPenDecreasing;
		for (int i = 8; i >= 0; i--){
			double cutoff = double(i)/10;
			//gapPenDecreasing = -3.-4.*i/10;
			multipleAlignment2ndRound=rawSequences.performMSAnextRound(prf,fprf,gapPen,gapExt,verboseMode,weightsModeOn,cutoff,codonLength,identities);
			//multipleAlignment2ndRound=rawSequences.performMSAnextRound(prf,fprf,gapPenDecreasing,gapExt,verboseMode,weightsModeOn,cutoff,domainScore,phosphScore,motifScore);
		}
		//txtProc::writeAlignmentWithoutCodeToFile(multipleAlignment2ndRound,encSeq,outputPrefix);						//write multiple alignment to a fileA
		//gapPenDecreasing = -2;
		multipleAlignment2ndRound=rawSequences.performMSAnextRound(prf,fprf,gapPen,gapExt,verboseMode,weightsModeOn,0,codonLength, identities);
		//multipleAlignment2ndRound=rawSequences.performMSAnextRound(prf,fprf,gapPenDecreasing,gapExt,verboseMode,weightsModeOn,0,codonLength, identities);
		txtProc::writeAlignmentToFile(multipleAlignment2ndRound,encSeq,outputPrefix);						//write multiple alignment to a fileA
	}
	else{
		std::cout << "input file and/or gap penalty not specified, try --help for more information" << std::endl;
	}
}
