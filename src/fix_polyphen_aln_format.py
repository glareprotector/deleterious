import global_stuff
import os
import pdb

all_files = os.listdir(global_stuff.polyphen_msa_directory)



for name in all_files:
    if name[-4:] == '.aln':
        print name
        f = open(global_stuff.polyphen_msa_directory + name, 'r')
        all_lines = f.readlines()
        g = open(global_stuff.polyphen_msa_directory + name + '.mine', 'w')
        #g.write(all_lines[0])
        #g.write(all_lines[1])
        for i in range(2, len(all_lines)):

            first = all_lines[i][0:67]
            new_first = first.split()[0]
            second = all_lines[i][70:]
            g.write('>' + first + '\n')
            g.write(second)
        f.close()
        g.close()
        
