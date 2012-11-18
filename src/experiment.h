
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
  string get_roc_file_name();

  ih m_options;
  dh m_parameters;

  experiment(string protein_list_file, string base_folder, string results_folder, ih options, dh parameters);
  string get_base_folder();
  int get_num_proteins();

  void get_roc_file();

};

#endif
