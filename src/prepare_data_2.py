import os
import string

import global_stuff

import sys
"""
deleterious mutation file
neutral mutation file
output protein list file
"""
protein_list_file = sys.argv[3]
deleterious_mutations_file = sys.argv[1]
neutral_mutations_file = sys.argv[2]

def get_mutations(mutations_file, which, the_dict):
    f = open(mutations_file)
    for line in f:
        s = line.strip().split('\t')
        protein_name = s[0]
        to_add = (s[1], s[2], s[3], which)
        #print protein_name, to_add
        try:
            the_dict[protein_name].append(to_add)
        except KeyError:
            the_dict[protein_name] = [to_add]
    return the_dict

def write_mutations(the_list, the_file):
    f = open(the_file, 'w')
    for elt in the_list:
        f.write(elt[0] + '\t' + elt[1] + '\t' + elt[2] + '\t' + str(elt[3]) + '\n')
    f.close()

the_dict = {}

get_mutations(deleterious_mutations_file, 0, the_dict)
import pdb
pdb.set_trace()
get_mutations(neutral_mutations_file, 1, the_dict)

protein_list = the_dict.keys()



f = open(global_stuff.all_seqs_file, 'r')

i = 0

for line in f:
    if i % 50 == 0:
        print i
    i += 1
    seq_name = string.split(line, sep = '|')[1]
    import pdb
    if seq_name not in protein_list:
        f.next()
    else:
        line0 = line
        line1 = f.next()

        folder = global_stuff.base_folder + seq_name + '/'
        try:
            os.makedirs(folder)
        except Exception, err:
            print err

        f1 = open(folder + 'seq', 'w')
        f1.write(line0)
        f1.write(line1)
        f1.close()

        mutations_file = folder + 'all_mutations'
        try:
            write_mutations(the_dict[seq_name], mutations_file)
        except KeyError:
            f2 = open(mutations_file, 'w')
            f2.close()

f2 = open(protein_list_file, 'w')

for name in protein_list:
    f2.write(name + '\n')

f2.close()
