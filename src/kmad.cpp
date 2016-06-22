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
    int refine_limit = 0;
    double ptm_modifier = 0;
    double strct_modifier = 0;
    double domain_modifier = 0;
    double motif_modifier = 0;
    double gap_ext_pen = 0;
    double gap_open_pen = 0;
    double end_pen = 0;
    bool out_encoded = false;
    bool one_round = false;
    bool first_gapped = false;
    bool refine = false;
    bool fade_out = false;
    bool optimize = false;
    bool no_feat = false;
    std::string filename;
    std::string output_prefix;
    std::string conf_file;
    std::string sbst_mat;
    double conservation_cutoff = 0;

    // TODO: Try to improve readability
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
       po::value<double>(&ptm_modifier)->default_value(10),
       "score for aligning phosphorylated residues")
      ("strct,s",
       po::value<double>(&strct_modifier)->default_value(4),
       "score for aligning secondary structure elements")
      ("d,d",
       po::value<double>(&domain_modifier)->default_value(3),
       "score for aligning domains")
      ("m,m",
       po::value<double>(&motif_modifier)->default_value(3),
                                 "probability multiplier for motifs")
      ("out-encoded",
       po::value<bool>(&out_encoded)->implicit_value(true)
                                    ->default_value(false),
       "output alignment with encoded features")
      ("end,n",
       po::value<double>(&end_pen)->default_value(-0.1),
       "penalty for gaps at the end (and beginning)")
      ("one-round, r", po::value<bool>(&one_round)->implicit_value(true)
                                                  ->default_value(false),
       "perform only one round of alignments (all against first)")
      ("feat_cutoff", po::value<double>(&conservation_cutoff)
          ->default_value(0.5),
       "conservation cutoff for the feature consensus")
      ("mat", po::value<std::string>(&sbst_mat)->default_value("DISORDER"),
       "substitution matrix")
      ("gapped", po::value<bool>(&first_gapped)->default_value(false)
                                               ->implicit_value(true),
       "'first sequence with gaps' mode"
       )
      ("no-feat", po::value<bool>(&no_feat)->default_value(false)
                                           ->implicit_value(true),
       "align without features"
       )
      ("refine", po::value<bool>(&refine)->default_value(false)
                                         ->implicit_value(true),
       "take alignment as input and refine it"
      )
      ("refine_limit", po::value<int>(&refine_limit)->default_value(0)
       )
      ("opt", po::value<bool>(&optimize)->default_value(false)
                                        ->implicit_value(true),
       "run alignment with optimizer"
       )
      ("fade-out", po::value<bool>(&fade_out)->default_value(false)
                                             ->implicit_value(true),
       "fade out feature scores with sequence identity"
       )
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

    // check if the obligatory parameters are provided
    if (vm.count("input") == 0 ||
        vm.count("gap_penalty") == 0 ||
        vm.count("output") == 0) {
      std::cout << "input file and/or gap penalty not specified, try --help \
                    for more information" << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // check if parameters are sane
    // TODO: Can boost program options check these rule for us?
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
    } else if (ptm_modifier < 0 || domain_modifier < 0
        || motif_modifier < 0 || strct_modifier < 0) {
      std::cout << "Scores for aligning features cannot be negative"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // Load custom feature settings

    f_config::FeatureSettingsMap f_set;
    // if the '--conf' option is chosen parse the configuration file
    if (vm.count("conf") == 1) {
      f_set = f_config::ConfParser::parse_conf_file(conf_file);
    }

    // Load sequence data

    fasta::FastaData fasta_data;
    try {
      fasta_data = fasta::parse_fasta(filename, codon_length, refine,
          refine_limit);
    } catch(const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    bool gapped = false;

    // Combine sequence and feature settings
    //
    auto sequence_data_plain = seq_data::process_fasta_data(
        fasta_data, f_set, gapped
    );

    // Perform the alignment

    // TODO: Maybe we can use a design pattern here because the type of
    // alignment depends on whether or not we will refine.
    std::vector<fasta::SequenceList> alignment;
    if (!refine) {
      alignment = msa::run_msa(sequence_data_plain,
                               f_set, gap_open_pen,
                               gap_ext_pen, end_pen, domain_modifier,
                               motif_modifier, ptm_modifier, strct_modifier,
                               codon_length,
                               one_round, sbst_mat, first_gapped, optimize,
                               fade_out, no_feat);
    } else {
      bool gapped = true;
      seq_data::SequenceData sequence_data_alignment = seq_data::process_fasta_data(
          fasta_data, f_set, gapped);
      if (refine_limit == 0) {
        refine_limit = fasta_data.sequences.size();
      }

      if (!check_length(fasta_data.sequences, refine_limit)) {
        throw std::runtime_error("In the 'refine' mode all sequences should "
                                 "have the same length");
      }
      alignment = msa::refine_alignment(sequence_data_plain,
                                        sequence_data_alignment,
                                        f_set, gap_open_pen,
                                        gap_ext_pen, end_pen, domain_modifier,
                                        motif_modifier, ptm_modifier,
                                        strct_modifier,
                                        codon_length,
                                        one_round, sbst_mat, first_gapped,
                                        optimize, fade_out, refine_limit,
                                        no_feat);
    }

    // Write alignment to file
    // TODO: al_out_index is always 1. Also, what is it?
    int al_out_index = 1;
    if (first_gapped) {
      first_gapped = 0;
    }

    // TODO: Probably can use a design pattern here. Encoding depends on user
    // option. Strategy might fit.
    if (out_encoded) {
      outfile::write_encoded_alignment(alignment[al_out_index],
                                       sequence_data_plain,
                                       output_prefix);
    } else {
      outfile::write_decoded_alignment(alignment[al_out_index],
                                       sequence_data_plain,
                                       output_prefix);
    }
}
