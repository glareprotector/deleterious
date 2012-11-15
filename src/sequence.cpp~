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

double sequence::predict_position(int pos, char mutant_res){
  string column = m_msa.get_column(pos);
  int count = 0;
  int total = 0;
  for(int i = 0; i < m_msa.get_num_seq(); i++){
    if(column[i] == mutant_res){
      count++;
    }
    if(column[i] != '-'){
      total++;
    }
  }
  return ((double)count + 1.0) / ((double)total + 1.0);
}

vector<double> sequence::predict_all_positions(){
 
  vector<double> ans;
  for(int i = 0; i < m_mutations.size(); i++){
    ans.push_back(predict_position(m_mutations[i].pos, m_mutations[i].mutant_res));
  }
  return ans;
}
