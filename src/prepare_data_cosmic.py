# for each mutation, record multiplicity, and also for the gene it is found in, how many mutations it has

import global_stuff
import pdb

cosmic_file = global_stuff.data_folder + 'CosmicMutantExport_v62_291112.tsv'
mutation_file = global_stuff.data_folder + 'cosmic_mutations'
gene_file = global_stuff.data_folder + 'cosmic_genes'

class mutation(object):

    def __init__(self, site, wild, mut):
        self.site = site
        self.wild = wild
        self.mut = mut
        self.count = None
        self.gene_mult = None

    def __hash__(self):
        return self.site.__hash__() + self.wild.__hash__() + self.mut.__hash__()

    def set_count(self, count):
        self.count = count

    def set_gene_mult(self, mult):
        self.gene_mult = mult

    def __lt__(self, other):
        return self.site < other.site

    def __str__(self):
        return str(self.site) + ','  + self.wild + ',' + self.mut + ',' + str(self.count) + ',' + str(self.gene_mult)


class site(object):
    def __init__(self, gene, pos):
        self.gene = gene
        self.pos = pos

    def __hash__(self):
        return self.gene.__hash__() + self.pos.__hash__()

    def __lt__(self, other):
        return self.gene < other.gene

    def __str__(self):
        return self.gene + ',' + str(self.pos)

    def __repr__(self):
        return self.__str__()

f = open(cosmic_file, 'r')

f.next()

mutation_counts = {}
gene_muts = {}
site_counts = {}
i = 0
for line in f:
    i += 1
    if i > 50000:
        break
    if i%50 == 0:
        print i

    try:
        s = line.strip().split('\t')
        gene = s[0]
        raw = s[14]
        raw = raw.split('.')[1]
        length = len(raw)
        wild = raw[0]
        mut = raw[length-1]
        pos = raw[1:(length-1)]
        a_site = site(gene, pos)
        to_add = mutation(a_site, wild, mut)
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
                gene_muts[gene].add(pos)
            except:
                gene_muts[gene] = set([pos])
    except Exception, err:
        x=2
        print err
        x=3

temp = [[key,mutation_counts[key]] for key in mutation_counts]
sorted_temp = sorted(temp, key = lambda x: x[0].site.gene)

for it in sorted_temp:
    print it[0]
    try:
        it[0].set_count(site_counts[it[0].site])
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
