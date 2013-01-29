f = open('../p53_old_2010/K8.data')
g = open('../data/train.mut')

i = 0
import pdb
h = open('../data/lathrop_single_p53.txt','w')

for mutline, featureline in zip(g,f):
    if len(mutline.strip().split('_')) == 1:
        h.write(mutline.strip() + '\t' + featureline.strip().split(',')[5408] + '\n')

h.close()
f.close()
g.close()
    

         
