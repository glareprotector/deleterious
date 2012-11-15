#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include "globals.h"
#include "msa.h"

class experiment;

class sequence{

 public:

  experiment* p_experiment;

  sequence(string _name, string _raw_seq, msa _msa);
  sequence(string _name, experiment* _p_experiment);
  sequence();

  string get_info_folder(string _name);
  vector<mutation> get_mutations();

  double predict_position(int pos, char mutant_res, ih options, dh parameters);
  vector<double> predict_all_positions(ih options, dh parameters);

  vector<double> get_sequence_weights(ih options, dh parameters);

  vector<mutation> m_mutations;
  string m_name;
  string m_raw_seq;
  msa m_msa;

};

#endif
