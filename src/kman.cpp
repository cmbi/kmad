#include "f_config.h"
#include "fasta.h"
#include "msa.h"
#include "outfile.h"
#include "seq_data.h"

#include <boost/program_options.hpp>
#include <ctime>
#include <iostream>
#include <stdexcept>


namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    int codon_length = 0;
    int ptm_modifier = 0;
    int domain_modifier = 0;
    int motif_modifier = 0;
    double gap_ext_pen = 0;
    double gap_open_pen;
    double end_pen = 0;
    bool out_encoded = false;
    bool one_round = false;
    std::string filename;
    std::string output_prefix;
    std::string conf_file;

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("input,i",
       po::value<std::string>(&filename),
       "input file name")
      ("output,o",
       po::value<std::string>(&output_prefix),
       "output file prefix")
      ("gap_penalty,g",
       po::value<double>(&gap_open_pen)->default_value(-5),
       "gap opening penalty")
      ("gap_extension,e",
       po::value<double>(&gap_ext_pen)->default_value(-1.),
       "gap extension penalty")
      ("codon_length,c",
       po::value<int>(&codon_length)->implicit_value(7)
                                    ->default_value(1),"codon length")
      ("phosph,p",
       po::value<int>(&ptm_modifier)->default_value(0),
       "score for aligning phosphorylated residues")
      ("domain,d",
       po::value<int>(&domain_modifier)->default_value(0),
       "score for aligning domains")
      ("motif,m",
       po::value<int>(&motif_modifier)->default_value(0),
                                 "probability multiplier for motifs")
      ("out-encoded",
       po::value<bool>(&out_encoded)->implicit_value(true)
                                    ->default_value(false),
       "Output alignment with encoded features")
      ("end,n",
       po::value<double>(&end_pen)->default_value(-0.1),
       "penalty for gaps at the end (and beginning)")
      ("one-round, r", po::value<bool>(&one_round)->implicit_value(true)
                                                  ->default_value(false),
       "perform only one round of alignments (all against first)")
      ("conf",
       po::value<std::string>(&conf_file),
       "configure file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << std::endl;
      std::exit(EXIT_SUCCESS);
    }

    if (vm.count("input") == 0 ||
        vm.count("gap_penalty") == 0 ||
        vm.count("output") == 0) {
      std::cout << "input file and/or gap penalty not specified, try --help \
                    for more information" << std::endl;
      std::exit(EXIT_FAILURE);
    }

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
    } else if (ptm_modifier < 0 || domain_modifier < 0 || motif_modifier < 0) {
      std::cout << "Scores for aligning features cannot be negative"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    
    f_config::FeatureSettingsMap f_set;
    if (vm.count("conf") == 1) {
      f_set = f_config::ConfParser::parse_conf_file(conf_file);
    }

    time_t start = clock();

    fasta::FastaData fasta_data;
    try {
      fasta_data = fasta::parse_fasta(filename, codon_length);
    } catch(const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }

    seq_data::SequenceData sequence_data = seq_data::process_fasta_data(
        fasta_data, f_set);

    std::cout << "runMsa" << std::endl;

    auto alignment = msa::run_msa(sequence_data, 
                                  f_set, gap_open_pen,
                                  gap_ext_pen, end_pen, domain_modifier, 
                                  motif_modifier, ptm_modifier, codon_length,
                                  one_round);

    if (out_encoded) {
      outfile::write_encoded_alignment(alignment[1], sequence_data,
                                       output_prefix);
    } else {
      outfile::write_decoded_alignment(alignment[1], sequence_data,
                                       output_prefix);
    }

    time_t end = clock();

    std::cout << "time: " << double(end - start)/CLOCKS_PER_SEC << std::endl;
}
