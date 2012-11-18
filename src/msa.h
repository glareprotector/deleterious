#ifndef MSA_H
#define MSA_H
#include <string>
#include <vector>

using namespace std;

class msa{
 public:
  vector<string> m_data;
  int length;
  int num_seq;

  msa(vector<string> _msa);
  msa();
  string get_column(int pos);
  string get_seq(int which);
  char operator()(int i, int j);

  int get_length();
  int get_num_seq();

  double get_distance(int pos1, int pos2);

};

ostream& operator<<(ostream& output, msa& the_msa);


struct mutation{
  int pos;
  char wild_res;
  char mutant_res;
  int is_deleterious;
  string name;
  mutation(int, char, char, int, string);
};

struct constraint{
  int pos;
  int res;
  bool satisfies(string seq);
  constraint(int _pos, char _res);
  constraint();
};

#endif
