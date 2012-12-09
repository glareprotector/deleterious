import global_stuff
import subprocess
import os

f = open(global_stuff.cosmic_gene_list, 'r')

for line in f:
    gene = line.strip()
    folder = global_stuff.COSMIC_BASE_FOLDER + gene + '/'
    try:
        os.makedirs(folder)
    except Exception, err:
        print err
    # find the sequence in cosmic downloaded seq folder
    seq_folder = global_stuff.cosmic_raw_data_folder + gene[0].upper() + '/'
    try:
        os.make_dirs(seq_folder)
    except Exception, err:
        pass
    seq_file = seq_folder + gene + '_protein.txt'
    new_seq_file = folder + 'seq'
    subprocess.call(['cp', seq_file, new_seq_file])

f.close()
