import sys
import os
import global_stuff

home = global_stuff.base_folder

input = sys.argv[1]

out = sys.argv[2]

f = open(input,'r')
g = open(out,'w')

i=0
for line in f:
    name = line.strip()
    print name, i
    i += 1
    if os.path.isdir(home + name):
        g.write(name + '\n')

f.close()
g.close()
