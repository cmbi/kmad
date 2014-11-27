#ifndef MSA_H
#define MSA_H

#include <iostream>
#include <vector>

class Sequences;
namespace msa{
  std::vector<std::string> run_msa(Sequences, std::string, double, double, 
                                   double, double, int, int, int, int, bool,
                                   std::vector<std::string>, std::vector<double>);
}

#endif /* MSA_H */
