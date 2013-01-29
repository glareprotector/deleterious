
import wc
import param
import objects
import global_stuff
import helper
import wrapper
import sys

name = sys.argv[1]
which_msa = int(sys.argv[2])
try:
    itera = int(sys.argv[3])
except:
    pass

p = param.param({'pdb':'1JOS', 'chain':'A', 'which_dataset':'CBS', 'uniprot_id':name, 'co':7.0, 'which_blast':0, 'which_msa':which_msa, 'ev':.05, 'blmax':999999,'hhblits_iter':itera, 'which_neighbors':1, 'protein_list_file':'rascalled_completed', 'to_leon':0, 'to_cluster':1, 'to_rascal':0, 'to_normd':0, 'norm_co':9.0, 'psiblast_iter':itera})

wc.get_stuff(wrapper.my_msa_obj_wrapper, p)

p.set_param('to_rascal', 1)

wc.get_stuff(wrapper.my_msa_obj_wrapper, p)

p.set_param('to_normd', 1)

wc.get_stuff(wrapper.my_msa_obj_wrapper, p)

