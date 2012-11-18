#include <string>
#include <vector>
#include "msa.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include "globals.h"
#include <cmath>

using namespace std;

msa read_msa(string file){

  ifstream f(file.c_str());
  string first_line;
  getline(f, first_line);
  istringstream s(first_line);
  string num_lines_string;
  getline(s, num_lines_string, '\t');
  int num_lines = atoi(num_lines_string.c_str());
  vector<string> ans;
  for(int i = 0; i < num_lines; i++){
    string temp;
    getline(f, temp);
    ans.push_back(temp);
  }

  return msa(ans);

}

string read_raw_seq(string file){

  ifstream f(file.c_str());
  string first_line;
  getline(f, first_line);
  string ans;
  getline(f, ans);

  return ans;
}

vector<mutation> read_mutation_file(string file, string name, int is_deleterious){

  vector<mutation> mutations;

  ifstream f(file.c_str());
  string temp;
  while(getline(f,temp)){
    istringstream stream(temp);
    string temp2;
    getline(stream, temp2, '\t');
    int pos = atoi(temp2.c_str()) - 1;
    getline(stream, temp2, '\t');
    char wild_res = temp2.c_str()[0];
    getline(stream, temp2);
    char mutant_res = temp2.c_str()[0];
    
    mutations.push_back(mutation(pos, wild_res, mutant_res, is_deleterious, name));
  }
  
  return mutations;
}

double string_distance(string x, string y){
  int num_diff = 0;
  int length = x.length();
  for(int i = 0; i < length; i++){
    if(x[i] != y[i]){
      num_diff++;
    }
  }
  return (double)num_diff / (double)length;
}
  

vector<string> read_protein_list_file(string file){

  vector<string> ans;
  ifstream f(file.c_str());
  while(f.good()){
    string temp;
    getline(f, temp);
    ans.push_back(temp);
  }
  return ans;
}

cdmap get_aa_weight_map(){
  cdmap ans;
  for(int i = 0; i < globals::all_aa.length(); i++){
    ans[globals::all_aa[i]] = 0.0;
  }
  return ans;
}

cimap get_aa_count_map(){
  cimap ans;
  for(int i = 0; i < globals::all_aa.length(); i++){
    ans[globals::all_aa[i]] = 0;
  }
  return ans;
}

bool pair_compare1(pair<double, int> x, pair<double, int> y){
  return x.first < y.first;
}

bool pair_compare2(pair<double, pair<int,int> > x, pair<double, pair<int,int> > y){
  return x.first < y.first;
}

bool satisfies_constraint_list_and(string seq, vector<constraint> constraint_list){
  for(int i = 0; i < constraint_list.size(); i++){
    if(!constraint_list[i].satisfies(seq)){
      return false;
    }
  }
  return true;
}

vector< vector<double> > read_double_mat(string file){
  string line;
  ifstream f(file);
  vector< vector<double>> ans;
  while(getline(f, line)){
    stringstream ss(line);
    string bit;
    vector<double> temp;
    while(getline(ss, bit, ',')){
      temp.push_back(atof(bit.c_str()));
    }
    ans.push_back(temp);
  }
  return ans;
}

double kl_distance(string x, string y){
  cdmap x_count;
  cdmap y_count;
  ccdmap xy_count;

  assert(x.length() == y.length());
  double len = x.length();
  for(int i = 0; i < x.length(); i++){
    x_count[x[i]] += 1.0 / len;
    y_count[y[i]] += 1.0 / len;
    xy_count[pair<char,char>(x[i],y[i])] += 1.0 / len;
  }
  double kl = 0;
  for(auto it = xy_count.begin(); it != xy_count.end(); ++it){
    //cout<<it->second<<" "<<x_count[it->first.first]<<" "<<y_count[it->first.second]<<endl;
    kl += it->second * log(it->second) - log(x_count[it->first.first]) - log(y_count[it->first.second]);
  }
  //cout<<kl<<endl;
  //exit(1);
  return kl;
}

bool display_string(string s){
  for(int i = 0; i < s.length(); i++){
    cout<<s[i]<<" ";
  }
  cout<<endl;
}
