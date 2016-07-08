#ifndef OUTFILE_H
#define OUTFILE_H

#include "fasta.h"


// TODO: opposite of what fasta does. One write actual fasta, other
//       writes custom fasta.

namespace outfile{
  void write_encoded_alignment(const fasta::SequenceList& sequences,
                               const fasta::FastaData& fasta_data,
                               const std::string& filename_prefix);
  void write_decoded_alignment(const fasta::SequenceList& sequences,
                               const fasta::FastaData& fasta_data,
                               const std::string& filename_prefix);
}

#endif /* OUTFILE_H */
