import global_stuff
import os
import pdb
import sys
import helper
"""
in_list absolute path
out_list absolute path
optional specifying params
"""

in_list = sys.argv[1]
out_list = sys.argv[2]

searching_for = sys.argv[3:]

f = open(in_list, 'r')

p = global_stuff.get_param()

helper.parse_p_input(p, sys.argv[3:])

completed = []

import objects

#to_check_for = [objects.general_distance, objects.general_msa, objects.neighbors_w_weight_w]
#to_check_for = [objects.neighbors_]
import wrapper
to_check_for = [wrapper.my_msa_obj_wrapper]

global_stuff.whether_to_look_at_whether_to_override = False

import wc
i = 0
for line in f:
    name = line.strip()
    if i % 50 == 0:
        print i
    i += 1

    p.set_param('uniprot_id',name)

    ok = True
    for obj in to_check_for:
        for filter_co in [0.35]:
            p.set_param('filter_co', filter_co)
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
