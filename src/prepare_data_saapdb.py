import global_stuff
import os
import urllib, urllib2
import string

f = open(global_stuff.saapdb_mutations_file,'r')

p_to_mut_list = {}

for line in f:
    s = line.strip().split(',')
    uniprot_id = s[2]
    wild_res = s[5]
    mut_res = s[31]
    pos = int(s[4])
    cls = int(s[78])
    mutation = [uniprot_id, pos, wild_res, mut_res, cls]
    try:
        p_to_mut_list[uniprot_id].append(mutation)
    except KeyError:
        p_to_mut_list[uniprot_id] = [mutation]

i = 0
bad = 0

h = open('saapdb_genes', 'w')
for u in p_to_mut_list:
    h.write(u + '\n')
h.close()




for uniprot_id in p_to_mut_list:
    directory = global_stuff.base_folder + uniprot_id + '/'
    try:
        os.makedirs(directory)
    except:
        pass

    print uniprot_id, i, len(p_to_mut_list), bad
    i += 1
    try:
        mutation_file = directory + 'all_mutations'
        g = open(mutation_file, 'w')
        for mutation in p_to_mut_list[uniprot_id]:
            g.write(string.join([str(m) for m in mutation],sep='\t') + '\n')
        g.close()

        url = 'http://www.uniprot.org/uniprot/' + uniprot_id + '.fasta'
        request = urllib2.Request(url,'')
        request.add_header('User-Agent', 'Python contact')
        response = urllib2.urlopen(request)
        page = response.read()

        seq_file = directory + 'seq'
        g = open(seq_file, 'w')
        g.write(page)
    except:
        bad += 1
    
"""
