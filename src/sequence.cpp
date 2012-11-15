#include "sequence.h"
#include <string>
#include "globals.h"
#include "helper.h"
#include "msa.h"
#include "experiment.h"
#include <iostream>

using namespace std;

sequence::sequence(string _name, string _raw_seq, msa _msa){
  m_name = _name;
  m_raw_seq = _raw_seq;
  m_msa = _msa;
}

string sequence::get_info_folder(string name){
  return p_experiment->get_base_folder() + name + string("/");
}

sequence::sequence(string _name, experiment* _p_experiment){

  p_experiment = _p_experiment;

  string info_folder = get_info_folder(_name);


  m_name = _name;
  m_raw_seq = read_raw_seq(info_folder + string("seq"));
  m_msa = read_msa(info_folder + string("msa"));

  vector<mutation> deleterious_mutations = read_mutation_file(info_folder + string("deleterious"), _name, 0);
  vector<mutation> neutral_mutations = read_mutation_file(info_folder + string("neutral"), _name, 1);
  
  m_mutations = deleterious_mutations;
  m_mutations.insert(m_mutations.end(), neutral_mutations.begin(), neutral_mutations.end());

}

sequence::sequence(){
}

vector<mutation> sequence::get_mutations(){
  return m_mutations;
}

double sequence::predict_position(int pos, char mutant_res, ih options, dh parameters){

  // most naive method.  can weigh sequences differently
  if(options[string("which_method")] == 0){
    vector<double> weights = get_sequence_weights(options, parameters);
    string column = m_msa.get_column(pos);
    double count = 0;
    double total = 0;
    for(int i = 0; i < m_msa.get_num_seq(); i++){
      if(column[i] == mutant_res){
	count += weights[i];
      }
      if(column[i] != '-'){
	total += weights[i];
      }
    }
  return (count + 1.0) / (total + 1.0);
  }
}

vector<double> sequence::get_sequence_weights(ih options, dh parameters){

  if(options[string("which_weight")] == 0){

    int num_seq_in_msa = m_msa.get_num_seq();
    vector< vector< double> > dists(num_seq_in_msa, vector<double>(num_seq_in_msa, 0));
    for(int i = 0; i < num_seq_in_msa; i++){
      for(int j = 0; j < num_seq_in_msa; j++){
	dists[i][j] = (m_msa.get_seq(i), m_msa.get_seq(j));
      }
    }
    vector<double> weights(num_seq_in_msa, 0);
    for(int i = 0; i < num_seq_in_msa; i++){
      double total;
      for(int j = 0; j < num_seq_in_msa; j++){
	total += dists[i][j];
      }
      weights[i] = (double)total / (double) num_seq_in_msa;
    }
    return weights;
  }
}

vector<double> sequence::predict_all_positions(ih options, dh parameters){
 
  vector<double> ans;
  for(int i = 0; i < m_mutations.size(); i++){
    ans.push_back(predict_position(m_mutations[i].pos, m_mutations[i].mutant_res, options, parameters));
  }
  return ans;
}
