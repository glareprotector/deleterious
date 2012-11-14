#include <string>
#include "experiment.h"
#include "helper.h"
#include <fstream>
#include <iostream>

using namespace std;


experiment::experiment(string protein_list_file, string _base_folder, string _results_folder):
  m_base_folder(_base_folder), m_results_folder(_results_folder){
  m_protein_names = read_protein_list_file(protein_list_file);
  m_proteins = vector<sequence>(0);
  for(int i = 0; i < m_protein_names.size(); i++){
    m_proteins.push_back(sequence(m_protein_names[i], this));
  }
}

void experiment::get_roc_file(){
  string roc_file = m_results_folder + string("roc");
  vector<double> all_scores;
  vector<mutation> all_mutations;
  for(int i = 0; i < m_proteins.size(); i++){
    vector<mutation> temp_mutations = m_proteins[i].get_mutations();
    all_mutations.insert(all_mutations.end(), temp_mutations.begin(), temp_mutations.end());
    vector<double> temp_scores = m_proteins[i].predict_all_positions();
    all_scores.insert(all_scores.end(), temp_scores.begin(), temp_scores.end());
  }

  ofstream f(roc_file.c_str());
  for(int i = 0; i < all_mutations.size(); i++){
    mutation temp = all_mutations[i];
    f<<all_scores[i]<<'\t'<<temp.is_deleterious<<'\t'<<temp.name<<'\t'<<temp.pos<<'\t'<<temp.wild_res<<'\t'<<temp.mutant_res<<'\n';
  }
  f.close();
}
  
string experiment::get_base_folder(){
  return m_base_folder;
}

int experiment::get_num_proteins(){
  return m_proteins.size();
}


