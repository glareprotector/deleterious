from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

import pdb

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
import helper
helper.parse_p_input(params, sys.argv[8:])




l = wc.get_stuff(objects.filtered_mutation_list_given_protein_list, params)


i = 0
my_l = []
for m in l:
    if i % size == rank:
        my_l.append(m)
    i += 1

import objects

print rank, len(my_l)

my_output = objects.get_output_obj(params, my_l, use_neighbor, ignore_pos, max_neighbor, num_trials, pseudo_total, helper.vanilla_similarity, helper.normalize_nothing, helper.mutation_to_class)

comm.Barrier()

root_output = comm.gather(my_output, root=0)

to_use = []

print rank, root_output
if rank == 0:
    for output in root_output:
        to_use = to_use + output

        helper.write_mat(to_use, output_file)