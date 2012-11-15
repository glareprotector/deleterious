#include <string>
#include <vector>
#include "msa.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>

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
    int pos = atoi(temp2.c_str());
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
