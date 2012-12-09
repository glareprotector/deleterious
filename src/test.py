import wc
import param
import objects
p = param.param({'ev':1e-10, 'protein_list_file':'hum_var_msa_dist_completed', 'uniprot_id':'P80075', 'avg_deg':20, 'n_cutoff':0, 'f_cutoff':15})
m = wc.get_stuff(objects.pairwise_dist, p, False, False, False)
