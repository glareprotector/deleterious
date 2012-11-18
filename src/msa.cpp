#include <string>
#include "msa.h"
#include <vector>
#include <iostream>

using namespace std;

msa::msa(vector<string> _msa){
  m_data = _msa;
}

string msa::get_column(int pos){

  char* ans = new char[get_num_seq()];

  for(int i = 0; i < get_num_seq(); i++){
    ans[i] = m_data[i][pos];
  }
  return string(ans);
}

int msa::get_num_seq(){
  return m_data.size();
}

int msa::get_length(){
  return m_data[0].size();
}

char msa::operator()(int i, int j){
  return m_data[i][j];
}

string msa::get_seq(int which){
  return m_data[which];
}



ostream& operator<<(ostream& output, msa& the_msa){
  for(int i = 0; i < the_msa.get_num_seq(); i++){
    for(int j = 0; j < the_msa.get_length(); j++){
      output<<the_msa(i,j);
    }
    output<<endl;
  }
  return output;
}

msa::msa(){
  m_data = vector<string>(0);
}

mutation::mutation(int _pos, char _wild_res, char _mutant_res, int _is_deleterious, string _name):
  pos(_pos), wild_res(_wild_res), mutant_res(_mutant_res), is_deleterious(_is_deleterious), name(_name){}

bool constraint::satisfies(string seq){
  return seq[pos] == res;
}

constraint::constraint(int _pos, char _res)
  :pos(_pos), res(_res){}

constraint::constraint(){
}
