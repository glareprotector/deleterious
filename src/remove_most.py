import os
import pdb

base_folder = '/home/fw27/d/deleterious/data/proteins/humvar/'

to_keep = ['seq', 'all_mutations']
i = 0
for a_dir in os.listdir(base_folder):
    print a_dir, i
    i += 1
    try:
        for a_file in os.listdir(base_folder+a_dir):
            full_path = base_folder + a_dir

            if a_file not in to_keep:
                os.remove(full_path + '/' + a_file)
    except Exception, err:
        print err
