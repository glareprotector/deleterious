
#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <string>
#include <vector>
#include "sequence.h"
#include <ext/hash_map>


using namespace std;

class experiment{

 public:

  vector<string> m_protein_names;
  vector<sequence> m_proteins;
  string m_base_folder;
  string m_results_folder;

  hash_map<string, int> options;
  hash_map<string, double> parameters;

  experiment(string protein_list_file, string base_folder, string results_folder);
  string get_base_folder();
  int get_num_proteins();

  void get_roc_file();

};

#endif
