#include "features_profile.h"
#include "profile.h"
#include "residue.h"
#include "scoring_matrix.h"
#include "sequences.h"
#include "misc.h"
#include "msa.h"
#include "txtproc.h"
#include "vec_util.h"

#include <boost/program_options.hpp>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <string>
namespace po = boost::program_options;
int main(int argc, char *argv[]){
  	int codon_length, phosphScore, domainScore, motifScore = 0;
  	double gapExt, gapPen, endPenalty, lcr_mod = 0;
  	bool out_encoded = false;
  	std::string filename, outputPrefix, conf_file;
  	po::options_description desc("Allowed options");
  	desc.add_options()
  		("help,h", "produce help message")
  		("input,i", po::value<std::string>(&filename), "input file name")
  		("output,o", po::value<std::string>(&outputPrefix), 
                                          "output file prefix")
  		("gap_penalty,g", po::value<double>(&gapPen)->default_value(-5),
                                          "gap opening penalty")
  		("gap_extension,e",po::value<double>(&gapExt)->default_value(-1.),
                                           "gap extension penalty")
  		("codon_length,c", po::value<int>(&codon_length)->implicit_value(7)
                                        ->default_value(1),"codon length")
  		("phosph,p", po::value<int>(&phosphScore)->default_value(0),
                                  "score for aligning phosphorylated residues")
  		("domain,d", po::value<int>(&domainScore)->default_value(0),
                                  "score for aligning domains")
  		("lcr,l", po::value<double>(&lcr_mod)->default_value(1),
                                  "gap penalty modifier inside a low complexity region")
  		("motif,m", po::value<int>(&motifScore)->default_value(0),
                                 "probability multiplier for motifs")
      ("out-encoded", po::value<bool>(&out_encoded)->implicit_value(true)
                                      ->default_value(false),
                                      "Output alignment with encoded features")
  		("end,n", po::value<double>(&endPenalty)->default_value(-0.1),
                                  "penalty for gaps at the end (and beginning)")
 		  ("conf", po::value<std::string>(&conf_file), "configure file");
  	po::variables_map vm;
  	po::store(po::parse_command_line(argc, argv, desc), vm);
  	po::notify(vm);
  	if (vm.count("help")) {
      		std::cout << desc << std::endl;
      		return 1;
  	}
  	if (vm.count("input") && vm.count("gap_penalty") 
        && vm.count("output")
        && misc::checkParameters(codon_length, phosphScore, domainScore,
                                 motifScore, gapExt, gapPen, endPenalty)){
  		time_t start = clock();
  		IDsList motifs_ids;
  		ProbsList motifs_probs;
      Sequences sequences;
      try{
          sequences = txtProc::read_fasta(filename, 
                                      codon_length, 
                                      &motifs_ids, 
                                      &motifs_probs);
      }
      catch(const std::exception& e)
      {
        std::cout << "Exception: " << e.what() << "\n";
        return -1;
      }
      StringSequences seq_names = sequences.get_names();
      StringSequences alignment = msa::run_msa(sequences, conf_file,
                                                gapPen, gapExt,
                                                endPenalty, lcr_mod, 
                                                domainScore,
                                                motifScore, 
                                                phosphScore,
                                                codon_length, 
                                                motifs_ids, 
                                                motifs_probs);
      if (out_encoded){
  		  txtProc::writeAlignmentToFile(alignment, seq_names, outputPrefix);						
      }
      else{
  		  txtProc::writeAlignmentWithoutCodeToFile(alignment, seq_names, 
                                                 outputPrefix, codon_length);						
      }
  		time_t end = clock();
  		std::cout << "time: " << double(end - start)/CLOCKS_PER_SEC << std::endl;
  	}
  	else{
  		std::cout << "input file and/or gap penalty not specified, try --help for more information" << std::endl;
  	}
}

