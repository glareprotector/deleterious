#include <string>
#include "experiment.h"
#include "helper.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>

using namespace std;


experiment::experiment(string protein_list_file, string _base_folder, string _results_folder, ih _options, dh _parameters):
  m_base_folder(_base_folder), m_results_folder(_results_folder){
  m_options = _options;
  m_parameters = _parameters;
  m_protein_names = read_protein_list_file(protein_list_file);
  m_proteins = vector<sequence>(0);
  for(int i = 0; i < m_protein_names.size(); i++){
    m_proteins.push_back(sequence(m_protein_names[i], this));
  }
}

string experiment::get_roc_file_name(){
  stringstream ss;
  ss << m_results_folder;
  for(auto it = m_options.begin(); it != m_options.end(); ++it){
    ss<<it->first + string("_")<<it->second<<"_";
  }
  for(auto it = m_parameters.begin(); it != m_parameters.end(); ++it){
    ss<<it->first + string("_")<<it->second<<"_";
  }
  return ss.str();
}

void experiment::get_roc_file(){
  string roc_base = get_roc_file_name();
  vector<double> all_scores;
  vector<mutation> all_mutations;
  for(int i = 0; i < m_proteins.size(); i++){
    cout<<"predicting: "<<m_proteins[i].m_name<<endl;
    vector<mutation> temp_mutations = m_proteins[i].get_mutations();
    all_mutations.insert(all_mutations.end(), temp_mutations.begin(), temp_mutations.end());
    vector<double> temp_scores = m_proteins[i].predict_all_positions(m_options, m_parameters);
    all_scores.insert(all_scores.end(), temp_scores.begin(), temp_scores.end());
  }
  string roc_file_name = roc_base + string(".roc");
  ofstream f(roc_file_name.c_str());
  
  for(int i = 0; i < all_mutations.size(); i++){
    mutation temp = all_mutations[i];
    f<<all_scores[i]<<'\t'<<temp.is_deleterious<<'\t'<<temp.name<<'\t'<<temp.pos<<'\t'<<temp.wild_res<<'\t'<<temp.mutant_res<<'\n';
  }
  f.close();
  string roc_plot_file_name = roc_base + string(".roc.pdf");
  string cmd = string("Rscript get_roc.r ") + roc_file_name + string(" ") + roc_plot_file_name;
  system(cmd.c_str());


}
  
string experiment::get_base_folder(){
  return m_base_folder;
}

int experiment::get_num_proteins(){
  return m_proteins.size();
}


