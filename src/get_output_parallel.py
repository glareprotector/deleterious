from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

import pdb

import sys
import wc, objects

output_file = sys.argv[1]
use_neighbor = sys.argv[2] == 'T'
ignore_pos = sys.argv[3] == 'T'
max_neighbor = int(sys.argv[4])
num_trials = int(sys.argv[5])
pseudo_total = float(sys.argv[6])
to_neighbor_p_value = sys.argv[7] == 'T'

import global_stuff
params = global_stuff.get_param()


import helper
helper.parse_p_input(params, sys.argv[8:])




#l = wc.get_stuff(objects.filtered_mutation_list_given_protein_list, params)
l = wc.get_stuff(objects.filtered_mutation_list, params)

i = 0
my_l = []
for m in l:
    if i % size == rank:
        my_l.append(m)
    i += 1

import objects

print rank, len(my_l)

which_dataset = params.get_param('which_dataset')

if which_dataset == 'cosmic' or which_dataset == 'their_cosmic':
    mutation_to_class = helper.temp_cosmic
elif which_dataset == 'humvar':
    mutation_to_class = helper.mutation_to_class
elif which_dataset == 'saapdb':
    mutation_to_class = helper.saapdb_to_class
elif which_dataset == 'p53':
    mutation_to_class = helper.p53_to_class
my_output = objects.get_output_obj(params, my_l, use_neighbor, ignore_pos, max_neighbor, num_trials, pseudo_total, helper.vanilla_similarity, helper.normalize_nothing, mutation_to_class, to_neighbor_p_value)

comm.Barrier()

root_output = comm.gather(my_output, root=0)

to_use = []


if rank == 0:
    for output in root_output:
        to_use = to_use + output

        helper.write_mat_raw(to_use, output_file)
