import os
import string

import global_stuff

def get_mutations(mutations_file):
    f = open(mutations_file)
    protein_to_mutations_dict = {}
    the_dict = {}
    for line in f:
        s = line.strip().split('\t')
        protein_name = s[0]
        to_add = (s[1], s[2], s[3])
        try:
            the_dict[protein_name].append(to_add)
        except KeyError:
            the_dict[protein_name] = [to_add]
    return the_dict

def write_mutations(the_list, the_file):
    f = open(the_file, 'w')
    for elt in the_list:
        f.write(elt[0] + '\t' + elt[1] + '\t' + elt[2] + '\n')
    f.close()

deleterious_mutations = get_mutations(global_stuff.deleterious_mutations_file)
neutral_mutations = get_mutations(global_stuff.neutral_mutations_file)

protein_list = set(deleterious_mutations.keys()) | set(neutral_mutations.keys())



f = open(global_stuff.all_seqs_file, 'r')

i = 0

for line in f:
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

        deleterious_file = folder + 'deleterious'
        try:
            write_mutations(deleterious_mutations[seq_name], deleterious_file)
        except KeyError:
            f2 = open(deleterious_file, 'w')
            f2.close()

        neutral_file = folder + 'neutral'
        try:
            write_mutations(neutral_mutations[seq_name], neutral_file)
        except KeyError:
            f2 = open(neutral_file, 'w')
            f2.close()

        deleterious_file = folder + 'deleterious'
        try:
            write_mutations(neutral_mutations[seq_name], neutral_file)
        except KeyError:
            f2 = open(neutral_file, 'w')
            f2.close()


f2 = open(global_stuff.protein_list_file, 'w')

for name in protein_list:
    f2.write(name + '\n')

f2.close()
