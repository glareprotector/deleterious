import global_stuff
import subprocess
import os
import sys
"""
cosmic_gene_list
cosmic_base_folder - remember to add slash at the end
cosmic_raw_data_folder - remember to add slash at the end
"""

cosmic_gene_list = sys.argv[1]
base_folder = global_stuff.base_folder
cosmic_raw_data_folder = global_stuff.cosmic_raw_data_folder

f = open(cosmic_gene_list, 'r')

for line in f:
    gene = line.strip()
    folder = base_folder + gene + '/'
    try:
        os.makedirs(folder)
    except Exception, err:
        print err
    # find the sequence in cosmic downloaded seq folder
    seq_folder = cosmic_raw_data_folder + gene[0].upper() + '/'
    try:
        os.make_dirs(seq_folder)
    except Exception, err:
        pass
    seq_file = seq_folder + gene + '_protein.txt'
    new_seq_file = folder + 'seq'
    subprocess.call(['cp', seq_file, new_seq_file])

f.close()
