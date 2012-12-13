import global_stuff
import os
import pdb
import sys

"""
in_list absolute path
out_list absolute path
"""

in_list = sys.argv[1]
out_list = sys.argv[2]

searching_for = sys.argv[3:]

f = open(in_list, 'r')

p = global_stuff.get_param()

completed = []

import objects

#to_check_for = [objects.general_distance, objects.general_msa, objects.neighbors_w_weight_w]
to_check_for = [objects.general_distance]

import wc

for line in f:
    name = line.strip()

    p.set_param('uniprot_id',name)

    ok = True
    for obj in to_check_for:

        inst = wc.get_wrapper_instance(obj)

        if not inst.has(p, False):
            ok = False
            

    if ok:
       completed.append(name)
    

g = open(out_list, 'w')
for name in completed:
    g.write(name + '\n')

f.close
g.close()
