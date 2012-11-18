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

  vector<int> get_neighbors(int pos, ih options, dh parameters);
  unordered_map< pair<int,int>, int> get_edge_to_rank(vector<vector<double>> dists);
  double get_distance(int pos1, int pos2);


  int get_length();

  vector<mutation> m_mutations;
  string m_name;
  string m_raw_seq;
  msa m_msa;
  vector< vector<double>> m_dists;
  unordered_map< pair<int, int>, int> m_edge_to_rank;
};



#endif
