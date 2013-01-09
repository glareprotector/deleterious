

import sys
import wc, objects
input_file = sys.argv[1]
output_file = sys.argv[2]
use_neighbor = sys.argv[3] == 'T'
ignore_pos = sys.argv[4] == 'T'
max_neighbor = int(sys.argv[5])
num_trials = int(sys.argv[6])
pseudo_total = float(sys.argv[7])

import global_stuff
params = global_stuff.get_param()

params.set_param('protein_list_file', input_file)
l = wc.get_stuff(objects.filtered_mutation_list_given_protein_list, params)




import objects
import helper
my_output = objects.get_output_obj(params, l, use_neighbor, ignore_pos, max_neighbor, num_trials, pseudo_total, helper.vanilla_similarity, helper.normalize_nothing, helper.mutation_to_class)





helper.write_mat(my_output, output_file)
