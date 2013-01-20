# for each mutation, record multiplicity, and also for the gene it is found in, how many mutations it has

import global_stuff
import pdb

#cosmic_file = global_stuff.data_folder + 'CosmicMutantExport_v62_291112.tsv'
#mutation_file = global_stuff.data_folder + 'cosmic_mutations_real'
#gene_file = global_stuff.data_folder + 'cosmic_genes_real'

cosmic_file = 'nras_raw'
mutation_file = 'NRAS_mutations'
gene_file = 'NRAS_genes'


f = open(cosmic_file, 'r')

f.next()

mutation_counts = {}
gene_muts = {}
site_counts = {}
i = 0
for line in f:
    i += 1
    if i%50 == 0:
        print i


    try:
        s = line.strip().split('\t')
        if s[15] == 'Substitution - Missense':

            gene = s[0]
            raw = s[14]
            raw = raw.split('.')[1]
            length = len(raw)
            wild = raw[0]
            mut = raw[length-1]
            pos = raw[1:(length-1)]
            a_site = helper.site(gene, pos)
            to_add = helper.mutation(a_site, wild, mut)

            if mut != wild and mut != '*':
                try:
                    mutation_counts[to_add] += 1
                except:
                    mutation_counts[to_add] = 1
                try:
                    site_counts[a_site] += 1
                except:
                    site_counts[a_site] = 1
                try:
                    gene_muts[gene].add(to_add)
                except:
                    gene_muts[gene] = set([to_add])
    except Exception, err:
        x=2
        print err
        x=3



temp = [[key,mutation_counts[key]] for key in mutation_counts]
sorted_temp = sorted(temp, key = lambda x: x[0].site.gene)

for it in sorted_temp:
    
    try:
        it[0].set_count(it[1])
        it[0].set_site_count(site_counts[it[0].site])
        it[0].set_gene_mult(len(gene_muts[it[0].site.gene]))
    except:

        x=2

pdb.set_trace()

mutation_list = [x[0] for x in sorted_temp]

temp = gene_muts.keys()
sorted_temp = sorted(temp)
gene_list = sorted_temp

g = open(mutation_file, 'w')
for mut in mutation_list:
    g.write(str(mut) + '\n')
g.close()

g = open(gene_file, 'w')
for gene in gene_list:
    g.write(gene + '\n')

import subprocess, os

for gene in gene_muts:
    folder = global_stuff.base_folder + gene + '/'
    try:
        os.makedirs(folder)
    except Exception, err:
        print err
    seq_folder = global_stuff.cosmic_raw_data_folder + gene[0].upper() + '/'
    seq_file = seq_folder + gene + '_protein.txt'
    new_seq_file = folder + 'seq'
    subprocess.call(['cp', seq_file, new_seq_file])
    mutation_file = folder + 'all_mutations'
    h = open(mutation_file, 'w')
    for mutation in gene_muts[gene]:
        h.write(str(mutation) + '\n')
    h.close()
