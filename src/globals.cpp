#include <string>
#include "globals.h"

string globals::base_folder = string("../data/proteins/humdiv/");
string globals::protein_list_file = string("../data/completed_list");
string globals::results_folder = string("/home/fultonw/Dropbox/deleterious/");
string globals::all_aa = string("ARNDCEQGHILKMFPSTWYV-");
string globals::all_aa_no_skip = string("ARNDCEQGHILKMFPSTWYV");
cimap globals::all_aa_to_i;

void globals::init_globals(){
  for(int i = 0; i < globals::all_aa.length(); i++){
    all_aa_to_i[globals::all_aa[i]] = i;
  }
}

