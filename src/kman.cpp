#include "fasta.h"
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
int main(int argc, char *argv[]) {
    int codon_length, phosph_score, domain_score, motif_score = 0;
    double gap_ext_pen, gap_open_pen, end_pen = 0;
    bool out_encoded = false;
    std::string filename, output_prefix, conf_file;
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("input,i", po::value<std::string>(&filename), "input file name")
      ("output,o", po::value<std::string>(&output_prefix), 
                                          "output file prefix")
      ("gap_penalty,g", po::value<double>(&gap_open_pen)->default_value(-5),
                                          "gap opening penalty")
      ("gap_extension,e", po::value<double>(&gap_ext_pen)->default_value(-1.),
                                           "gap extension penalty")
      ("codon_length,c", po::value<int>(&codon_length)->implicit_value(7)
                                        ->default_value(1),"codon length")
      ("phosph,p", po::value<int>(&phosph_score)->default_value(0),
                                  "score for aligning phosphorylated residues")
      ("domain,d", po::value<int>(&domain_score)->default_value(0),
                                  "score for aligning domains")
      ("motif,m", po::value<int>(&motif_score)->default_value(0),
                                 "probability multiplier for motifs")
      ("out-encoded", po::value<bool>(&out_encoded)->implicit_value(true)
                                      ->default_value(false),
                                      "Output alignment with encoded features")
      ("end,n", po::value<double>(&end_pen)->default_value(-0.1),
                                  "penalty for gaps at the end (and beginning)")
      ("conf", po::value<std::string>(&conf_file), "configure file");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) {
          std::cout << desc << std::endl;
          return 1;
    }
    if (vm.count("input") && vm.count("gap_penalty") && vm.count("output")) {
      time_t start = clock();
      // Parameter check
      if (codon_length < 1 || codon_length > 10) {
        std::cout << "Codon length ('-c' flag) should be between 1 and 10"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      } else if (gap_open_pen >= 0 || gap_ext_pen >= 0) {
        std::cout << "Gap opening penalty (-g) and gap extension penalty (-e) \
                      need to be lower than 0" << std::endl;
        std::exit(EXIT_FAILURE);
      } else if (end_pen > 0) {
        std::cout << "End gap penalty value (-n) cannot be a positive number"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      } else if (phosph_score < 0 || domain_score < 0 || motif_score < 0) {
        std::cout << "Scores for aligning features cannot be negative"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      //
      IDsList motifs_ids;
      ProbsList motifs_probs;
      Sequences sequences;
      try {
          sequences = fasta::parse_fasta(filename, 
                                         codon_length, 
                                         &motifs_ids, 
                                         &motifs_probs);
      } catch(const std::exception& e) {
        std::cout << "Exception: " << e.what() << "\n";
        std::exit(EXIT_FAILURE);
      }
      StringSequences seq_names = sequences.get_names();
      StringSequences alignment = msa::run_msa(sequences, conf_file,
                                               gap_open_pen, gap_ext_pen,
                                               end_pen,
                                               domain_score,
                                               motif_score, 
                                               phosph_score,
                                               codon_length, 
                                               motifs_ids, 
                                               motifs_probs);
      if (out_encoded) {
        txtproc::writeAlignmentToFile(alignment, seq_names, output_prefix);            
      } else {
        txtproc::writeAlignmentWithoutCodeToFile(alignment, seq_names, 
                                                 output_prefix, codon_length);            
      }
      time_t end = clock();
      std::cout << "time: " << double(end - start)/CLOCKS_PER_SEC << std::endl;
    } else{
      std::cout << "input file and/or gap penalty not specified, try --help \
                    for more information" << std::endl;
    }
}

