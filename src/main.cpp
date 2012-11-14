#include "experiment.h"
#include <string>
#include <iostream>

using namespace std;

int main(){
  experiment test = experiment(globals::protein_list_file, globals::base_folder, globals::results_folder);
  test.get_roc_file();
}
