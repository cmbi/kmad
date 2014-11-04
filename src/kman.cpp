#include "Residue.h"
#include "Sequences.h"
#include "FeaturesProfile.h"
#include "Profile.h"
#include "txtProc.h"
#include "vecUtil.h"
#include "misc.h"
#include <boost/program_options.hpp>
#include <ctime>
#include <iostream>
#include <string>
namespace po = boost::program_options;
int main(int argc, char *argv[]){
	int codonLength, phosphScore,domainScore, motifScore;
	double gapExt, gapPen, endPenalty, lcr_mod;
	bool weightsModeOn;
	std::string filename,verboseMode,outputPrefix,conf_file;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "produce help message")
		("input,i", po::value<std::string>(&filename), "input file name")
		("output,o", po::value<std::string>(&outputPrefix), "output file prefix")
		("gap_penalty,g", po::value<double>(&gapPen), "gap opening penalty")
		("gap_extension,e",po::value<double>(&gapExt)->default_value(-1.),"gap extension penalty")
		("codon_length,c", po::value<int>(&codonLength)->implicit_value(7)->default_value(1),"codon length")
		("phosph,p", po::value<int>(&phosphScore)->default_value(0),"score for aligning phosphoryated residues")
		("domain,d", po::value<int>(&domainScore)->default_value(0),"score for aligning domains")
		("lcr,l", po::value<double>(&lcr_mod)->default_value(1), "gap penalty modifier inside a low complexity region")
		("motif,m", po::value<int>(&motifScore)->default_value(0),"probability multiplier for motifs")
		("end,n", po::value<double>(&endPenalty)->default_value(-0.1),"penalty for gaps at the end (and beginning)")
		("weights,w", po::value<bool>(&weightsModeOn)->implicit_value(true)->default_value(false),"all sequences contribute to the profile with weights(=similarity)")
		("conf", po::value<std::string>(&conf_file), "configure file")
		("verbose,v",po::value<std::string>(&verboseMode)->implicit_value("1")->default_value("0"),"verbose mode");
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if (vm.count("help")) {
    		std::cout << desc << std::endl;
    		return 1;
	}
	if (vm.count("input") && vm.count("gap_penalty") && misc::checkParameters(codonLength,phosphScore,domainScore,motifScore,gapExt,gapPen,weightsModeOn,endPenalty)){
		time_t start = clock();
		std::vector<std::string> motifs_ids;
		std::vector<double> motifs_probs, identities;
		Sequences rawSequences(txtProc::processFASTA(filename,codonLength, &motifs_ids, &motifs_probs));
		Profile prf;
		FeaturesProfile fprf(domainScore,phosphScore,motifScore,lcr_mod,motifs_ids,motifs_probs);
		if (!conf_file.empty()){
			txtProc::process_conf_file(conf_file,fprf,rawSequences);
		}
		//first round of the alignment - all vs 1st
		std::vector<std::string> multipleAlignment(rawSequences.performMSAfirstround(prf,fprf,gapPen,endPenalty,gapExt,verboseMode,weightsModeOn,codonLength,identities));
		std::vector<std::vector<std::vector<std::string> > > encSeq = rawSequences.getEncodedSequences();
		std::vector<std::string> alignment2ndRound;
		int prev_alignments = 0;
		for (int i = 8; i >= 0; i--){
			double cutoff = double(i)/10;
			//gapPen = -3.-4.*i/10; //decreasing gap penalties
			rawSequences.performMSAnextRounds(&alignment2ndRound,prf,fprf,gapPen,endPenalty,gapExt,verboseMode,weightsModeOn,cutoff,codonLength,identities, prev_alignments);
			//prev_alignments - number of alignments performed in the previous round - to omit this round if the number of aligned sequences is the same as in the previous round
		}
		//gapPen = -2; //decreasing gap penalties
		prev_alignments = 0;  // to align (again) all sequences to the profile
		rawSequences.performMSAnextRounds(&alignment2ndRound,prf,fprf,gapPen,endPenalty,gapExt,verboseMode,weightsModeOn,0,codonLength, identities, prev_alignments);
		txtProc::writeAlignmentToFile(alignment2ndRound,encSeq,outputPrefix);						//write multiple alignment to a fileA
		//txtProc::writeAlignmentWithoutCodeToFile(alignment2ndRound,encSeq,outputPrefix);						//write multiple alignment to a fileA
		time_t end = clock();
		std::cout << "time: " << double(end - start)/CLOCKS_PER_SEC << std::endl;
	}
	else{
		std::cout << "input file and/or gap penalty not specified, try --help for more information" << std::endl;
	}
}
