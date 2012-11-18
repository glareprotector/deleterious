#include "sequence.h"
#include <string>
#include "globals.h"
#include "helper.h"
#include "msa.h"
#include "experiment.h"
#include <iostream>
#include <algorithm>

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
  m_dists = read_double_mat(info_folder + string("dists"));

  m_edge_to_rank = get_edge_to_rank(m_dists);

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
    vector<int> neighbors = get_neighbors(pos, options, parameters);
    vector<constraint> constraint_list;
    for(int i = 0; i < neighbors.size(); i++){
      constraint_list.push_back(constraint(neighbors[i], m_raw_seq[neighbors[i]]));
    }
    //cout<<"mutant_res: "<<mutant_res<<" "<<"pos: "<<pos<<endl;
    //display_string(m_msa.get_column(pos));
    for(int i = 0; i < m_msa.get_num_seq(); i++){
      if(satisfies_constraint_list_and(m_msa.get_seq(i), constraint_list)){
	if(column[i] == mutant_res){
	  count += weights[i];
	}
	if(column[i] != '-'){
	  total += weights[i];
	}
      }
    }
    //cout<<count<<" "<<total<<endl;
    return (count + 1.0) / (total + 1.0);
  }
}

vector<double> sequence::get_sequence_weights(ih options, dh parameters){

  if(options[string("which_weight")] == 0){

    int num_seq_in_msa = m_msa.get_num_seq();
    vector< vector< double> > dists(num_seq_in_msa, vector<double>(num_seq_in_msa, 0));
    for(int i = 0; i < num_seq_in_msa; i++){
      for(int j = 0; j < num_seq_in_msa; j++){
	dists[i][j] = string_distance(m_msa.get_seq(i), m_msa.get_seq(j));
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

  else if(options[string("which_weight")] == 1){
    int num_seq_in_msa = m_msa.get_num_seq();
    vector<double> weights(num_seq_in_msa, 0.0);
    for(int i = 0; i < m_msa.get_length(); i++){
      cimap pos_counts = get_aa_count_map();
      string col = m_msa.get_column(i);
      for(int j = 0; j < m_msa.get_num_seq(); j++){
	pos_counts[col[j]] += 1;
      }
      for(int j = 0; j < m_msa.get_num_seq(); j++){
	weights[j] += 1.0 / pos_counts[col[j]];
      }
    }
    //for(int i = 0; i < weights.size(); i++) cout<<weights[i]<<" ";
    return weights;
  }

}

double sequence::get_distance(int pos1, int pos2){

  string x = m_msa.get_column(pos1);
  string y = m_msa.get_column(pos2);
  string x_no_skip, y_no_skip;
  for(int i = 0; i < m_msa.get_num_seq(); i++){
    if(m_msa(i, pos1) != '-' && m_msa(i, pos2) != '-'){
      x_no_skip.push_back(m_msa(i, pos1));
      y_no_skip.push_back(m_msa(i, pos2));
    }
  }
  return kl_distance(x_no_skip, y_no_skip);
}

unordered_map< pair<int,int>, int> sequence::get_edge_to_rank(vector< vector<double> > dists){
  unordered_map< pair<int,int>, int> edge_to_rank;
  vector< pair<double, pair<int,int> > > temp;
  for(int i = 0; i < dists.size(); i++){
    for(int j = 0; j < dists.size(); j++){
      temp.push_back(pair< double, pair<int,int>>(dists[i][j], pair<int,int>(i,j)));
    }
  }
  sort(temp.begin(), temp.end(), pair_compare2);
  int idx = temp.size() - 1;
  for(int i = 0; i < temp.size(); i++){
    edge_to_rank[temp[i].second] = idx;
    idx--;
  }
  return edge_to_rank;
}


vector<int> sequence::get_neighbors(int pos, ih options, dh parameters){
  if(options[string("which_neighbor_method")] == 0){
    int num_aa = globals::all_aa.length();
    vector<double> kls(get_length(), 0.0);
    //cout<<pos<<endl;
    for(int i = 0; i < get_length(); i++){
      vector< vector< double> > joints(num_aa, vector<double>(num_aa, 0.0));
      vector<double> x_marg(num_aa, 0.0);
      vector<double> y_marg(num_aa, 0.0);
      int count = 0;
      string x_no_skip, y_no_skip;
      for(int k = 0; k < m_msa.get_num_seq(); k++){
	if(m_msa(k, pos) != '-' && m_msa(k,i) != '-'){
	  x_no_skip.push_back(m_msa(k, pos));
	  y_no_skip.push_back(m_msa(k, i));
	}
      }
      if(x_no_skip.length() > m_msa.get_num_seq() / 2){
	kls[i] = kl_distance(x_no_skip, y_no_skip);
      }
      else{
	kls[i] = -1;
      }

    }
    vector< pair<double, int> > temp(get_length());
    for(int i = 0; i < get_length(); i++){
      temp[i] = pair<double, int>(kls[i], i);
    }
    sort(temp.begin(), temp.end(), pair_compare1);
    vector<int> ans(0);
    int i = 0;
    int idx = temp.size() - 1;
    while(i < parameters[string("num_neighbor")]){
      if(idx != pos){
	//	cout<<temp[idx].second<<" is neighbor with kl "<<temp[idx].first<<endl;
	//display_string(m_msa.get_column(temp[idx].second));
	ans.push_back(temp[idx].second);
	i++;
	idx--;
      }
      //      cout<<endl;
    }
    return ans;
  }
  else if(options[string("which_neighbor_method")] == 1){
    int x;
  }
}
    
    

vector<double> sequence::predict_all_positions(ih options, dh parameters){
 
  vector<double> ans;
  for(int i = 0; i < m_mutations.size(); i++){
    assert(m_mutations[i].wild_res == m_raw_seq[m_mutations[i].pos]);
    //ans.push_back(predict_position(m_mutations[i].pos, m_mutations[i].wild_res, options, parameters));
    ans.push_back(predict_position(m_mutations[i].pos, m_mutations[i].mutant_res, options, parameters));
  }
  return ans;
}

int sequence::get_length(){
  return m_raw_seq.length();
}
