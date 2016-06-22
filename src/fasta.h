#ifndef FASTA_H
#define FASTA_H
#include "f_config.h"

#include <iostream>
#include <unordered_map>
#include <vector>

/// Parses a fasta or semi-fasta(encoded) file to a vec of Sequence structs,
namespace fasta {
    /// Holds a single residue codon (x characters that encode the amino acid
    /// and it's features), and a list of its features
    struct Residue {
        Residue(std::string codon, std::vector<std::string> features)
            : codon(codon),
            features(features) {}
        explicit Residue(std::string codon) : codon(codon) {}
        Residue(){};

        std::string codon;
        std::vector<std::string> features;
    };

    /// Holds description which is later used as a fasta header and a list of
    /// residues
    struct Sequence {
        std::string description;
        std::vector<Residue> residues;
    };
    typedef std::vector<fasta::Sequence> SequenceList;

    /// Holds all the data from the (semi)fasta file - sequences with their
    /// headers, list of motif probabilities, list of all features present
    /// in the input
    struct FastaData {
        SequenceList sequences;
        std::unordered_map<std::string, double> probabilities;
        std::vector<std::string> feature_list;
    };

    /// \brief Parses a 'fasta' file
    ///
    /// \throws std::runtime_error when motif probabilities format is incorrect
    FastaData parse_fasta(std::string const& filename, int codon_length);


    /// Make a sequence
    Sequence make_sequence(const std::string& description,
                           const std::string& codons, int codon_length);

    /// Makes a residue
    Residue make_residue(const std::string& codon);

    /// Converts a Sequence struct to a string
    std::string sequence_to_string(const Sequence& seq);

    /// Check if all sequence lengths are equal (for the refinement mode)
    ///
    /// \param[in] limit The number of sequence to include in the length check.
    bool check_length(SequenceList const& sequences, int limit);

    FastaData get_conf_data(const FastaData& fasta_data,
                            const f_config::FeatureSettingsMap& f_set,
                            bool gapped);
    /// \brief Removes gaps from all sequences in the list
    SequenceList remove_gaps(const SequenceList& sequences);

    /// \brief add features based on the provided regular expression @pattern
    void assign_feature_by_pattern(fasta::SequenceList& sequences,
                                   const std::string& pattern,
                                   const std::string& feat_name);

    /// \brief find position in the alignment based on position in the sequence
    ///  e.g. for @sequence '--A' and @position 0 will return 2
    int find_real_pos(const std::string& sequence, int position);

    /// creates a list of all features present in @sequences + default ptm
    ///         features
    FeatureNamesList make_feature_list(const fasta::SequenceList& sequences);
}


#endif /* FASTA_H */
