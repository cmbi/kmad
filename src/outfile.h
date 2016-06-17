#ifndef OUTFILE_H
#define OUTFILE_H

#include "fasta.h"
#include "seq_data.h"


// TODO: opposite of what fasta does. One write actual fasta, other
//       writes custom fasta.

namespace outfile{
  void write_encoded_alignment(const fasta::SequenceList& sequences,
                               const seq_data::SequenceData& sequence_data,
                               const std::string& filename_prefix);
  void write_decoded_alignment(const fasta::SequenceList& sequences,
                               const seq_data::SequenceData& sequence_data,
                               const std::string& filename_prefix);
}

#endif /* OUTFILE_H */
