#include "experiment.h"
#include <string>
#include <iostream>
#include "globals.h"

using namespace std;

int main(){

  ih options;
  dh parameters;
  options[string("which_method")] = 0;
  options[string("which_weight")] = 1;
  options[string("which_neighbor_method")] = 1;
  for(int i = 1; i < 4; i++){
    parameters[string("edge_avg")] = 45*i;

    experiment test = experiment(globals::protein_list_file, globals::base_folder, globals::results_folder, options, parameters);
    test.get_roc_file();
  }
  /*for(int i = 0; i < 10; i++){
    parameters[string("num_neighbor")] = 1;
    experiment test = experiment(globals::protein_list_file, globals::base_folder, globals::results_folder, options, parameters);
    test.get_roc_file();
    }*/
}
