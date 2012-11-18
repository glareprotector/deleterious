#ifndef HELPER_H
#define HELPER_H

#include <string>
#include "msa.h"
#include <vector>
#include "globals.h"


msa read_msa(string file);
string read_raw_seq(string file);
vector<mutation> read_mutation_file(string file, string name, int is_deleterious);
vector<string> read_protein_list_file(string file);
vector<vector<double>> read_double_mat(string file);
double string_distance(string x, string y);
cimap get_aa_count_map();
cdmap get_aa_weight_map();
bool pair_compare1(pair<double, int> x, pair<double, int> y);
bool pair_compare2(pair<double, pair<int,int> > x, pair<double, pair<int,int> > y);
bool satisfies_constraint_list_and(string seq, vector<constraint> constraint_list);
void display_string(string s);
double kl_distance(string x, string y);



#endif
