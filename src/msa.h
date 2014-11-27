#ifndef MSA_H
#define MSA_H

#include "types.h"

#include <iostream>
#include <vector>

class Sequences;
namespace msa{
  string_sequences run_msa(Sequences, std::string, double, double, 
                           double, double, int, int, int, int, bool,
                           ids_list, probs_list);
}

#endif /* MSA_H */
