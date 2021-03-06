import string
#import Bio.PDB
import csv
#import constants
import string
import re
import math
import global_stuff
import pdb
import param
import random
import sys
from Bio.Align import MultipleSeqAlignment
import subprocess

class mutation(object):

    def __init__(self, site, wild, mut):
        self.site = site
        self.wild = wild
        self.mut = mut
        self.count = None
        self.gene_mult = None
        self.site_count = None

    def __hash__(self):
        return self.site.__hash__() + self.wild.__hash__() + self.mut.__hash__()

    def set_count(self, count):
        self.count = count

    def set_site_count(self, count):
        self.site_count = count

    def set_gene_mult(self, mult):
        self.gene_mult = mult

    def __lt__(self, other):
        return self.site < other.site

    def __eq__(self, other):
        return self.__str__() == other.__str__()

    def __str__(self):
        return str(self.site) + ','  + self.wild + ',' + self.mut + ',' + str(self.count) + ',' + str(self.site_count) + ',' + str(self.gene_mult)


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

    def __eq__(self, other):
        return self.__str__() == other.__str__()

    def __repr__(self):
        return self.__str__()


# also sets global_stuff.RESULTS_FOLDER to proper value
def read_param(file_location):
    
    # read folder_name
    f = open(constants.INFO_FOLDER + file_location)
    the_params = param.param({})
    hp_values = param.param()
    folder_name = f.readline().strip()

    global_stuff.RESULTS_FOLDER = global_stuff.RESULTS_BASE_FOLDER + folder_name + '/'

    for line in f.readlines():
        print >> sys.stderr, line
        if line[0] != '#':
            s = line.strip().split(',')
            if s[0] != 'hp':
                the_name = s[0]
                if the_name == 'n':
                    node_features = []
                    for i in range(1, len(s)):
                        node_features.append(constants.get_master_node_feature_list()[int(s[i])])
                    the_params.set_param('n', node_features)
                if the_name == 'e':
                    edge_features = []
                    for i in range(1, len(s)):
                        edge_features.append(constants.get_master_edge_feature_list()[int(s[i])])
                    the_params.set_param('e', edge_features)
                try:
                    the_type = s[1]
                    if the_type == 'f':
                        the_params.set_param(the_name, float(s[2]))
                    elif the_type == 'i':
                        the_params.set_param(the_name, int(s[2]))
                    elif the_type == 's':
                        the_params.set_param(the_name, s[2])
                except:
                    pass

    # hp values file happens to be the same as info file, so set that

    the_params.set_param('hpvf', file_location)


    if len(the_params.get_param('e')) != 0:
        assert the_params.get_param('wif') != 2

    
    return folder_name, the_params


# query should be a SeqRecord object
def filter_msa_based_on_pos_neighbors_and_query(query, pos, msa, neighbors):

    seqs = [query]
    s = query.seq.tostring()

    ans = ''

    for i in range(len(msa)):
        okay = True
        for n in neighbors:
            if msa[i][n] != s[n]:
                okay = False
                break
        if okay:
            ans += msa.get_column(pos)[i]
            seqs.append(msa[i])
            
    return MultipleSeqAlignment(seqs)

def filter_msa_based_on_pos_constraint(msa, constraints):

    seqs = []
    for i in range(len(msa)):
        works = True
        for c in constraints:
            pos = c[0]
            res = c[1]
            if msa[i][pos] != res:
                works = False
                break
        if works:
            seqs.append(msa[i])
    return MultipleSeqAlignment(seqs)


def get_weight_of_msa_seqs(msa):

    unnormalized = [0 for i in range(len(msa))]
    num = [1 for i in range(len(msa))]

    for i in range(msa.get_alignment_length()):
        d = {}
        col = msa.get_column(i)
        for j in range(len(msa)):
            """
            if col[j] != '-':
                num[j] += 1
                try:
                    d[col[j]] += 1
                except:
                    d[col[j]] = 1
            """
            num[j] += 1
            try:
                d[col[j]] += 1
            except:
                d[col[j]] = 1
        for j in range(len(msa)):
            """
            if col[j] != '-':
                unnormalized[j] += 1.0 / d[col[j]]
            """
            unnormalized[j] += 1.0 / d[col[j]]
    try:
        normalized = [unnormalized[i] / num[i] for i in range(len(msa))]
    except:
        x = 2
        import pdb
#        pdb.set_trace()
    return normalized

def predict_position(params, recalculate, mutation, use_neighbor):

    import wc
    import objects

    #use_neighbor = False

    protein_name = mutation[0]
    pos = mutation[1]
    params.set_param('uniprot_id', protein_name)
    seq = wc.get_stuff(objects.dW, params, recalculate, False, False)
    msa = wc.get_stuff(objects.agW, params, recalculate, False, False)

    if use_neighbor:
    
        neighbors = wc.get_stuff(objects.neighbors_w, params, recalculate, False, False)
        try:
            msa = filter_msa_based_on_pos_neighbors_and_query(seq, pos, msa, neighbors[pos])
        except:
            x = 2
            pdb.set_trace()
    weights = get_weight_of_msa_seqs(msa)

    res = mutation[3]

    score = 0

    col = msa.get_column(pos)
    for i in range(len(msa)):
        if col[i] == res:
            score += weights[i]

    return score / sum(weights)




def predict_position_energy(params, recalculate, mutation, use_neighbor, ignore_pos):

    import wc
    import objects

    #use_neighbor = False

    protein_name = mutation[0]
    pos = mutation[1]
    wild_res = mutation[2]
    mut_res = mutation[3]
    params.set_param('uniprot_id', protein_name)
    seq = wc.get_stuff(objects.dW, params, recalculate, False, False)
    msa = wc.get_stuff(objects.agW, params, recalculate, False, False)

    score = 0
    
    col = msa.get_column(pos)

    if not ignore_pos:
        try:
            score += -1.0 * math.log( float(col.count(mut_res)+1) / (col.count(wild_res)+1) )
            #score = -1.0 * ( float(col.count(mut_res)+1) / (col.count(wild_res)+1) )
        except:
            pdb.set_trace()
            x=2
        print >> sys.stderr, score
    if col.count(mut_res) == 0:
            return -1.0 * score
    
    
    if use_neighbor:
        constraints_a = [(pos,wild_res)]
        filter_a_msa = filter_msa_based_on_pos_constraint(msa, constraints_a)
        constraints_b = [(pos,mut_res)]
        filter_b_msa = filter_msa_based_on_pos_constraint(msa, constraints_b)
        all_neighbors = wc.get_stuff(objects.neighbors_w, params, recalculate, False, False)
        neighbors = all_neighbors[pos]
        for neighbor in neighbors[0:1]:
            try:
                na_col = filter_a_msa.get_column(neighbor)
                nb_col = filter_b_msa.get_column(neighbor)
            except:
                pdb.set_trace()
                x=2
            na_col_no_skip = [x for x in na_col if x != '-']
            nb_col_no_skip = [x for x in nb_col if x != '-']
            score += get_KL_real(nb_col_no_skip, na_col_no_skip)
            if abs(score) < .001:
                pdb.set_trace()
                x=2

    return -1.0 * score

def normalize(vect):
    total = 0
    for x in vect:
        total += x
    return [x/total for x in vect]

def vanilla_similarity(aa1, aa2):
    if aa1 not in global_stuff.ignore_aas and aa2 not in global_stuff.ignore_aas:
        if aa1 == aa2:
            return 1.0
        else:
            return 0.0
    else:
        return 0.0

def category_similarity(aa1, aa2):
    if aa1 not in global_stuff.ignore_aas and aa2 not in global_stuff.ignore_aas:
        if global_stuff.aa_to_class[aa1] == global_stuff.aa_to_class[aa2]:
            return 1.0
        else:
            return 0.0
    else:
        return 0.0

# distance function should return something between 0 and 1.  that means distance to each residue can have a max of 1
def compute_similarity_score_to_residue(column, weights, residue, similarity_f):

    similarity = 0
    for i in range(len(column)):
        
        similarity += weights[i] * similarity_f(residue, column[i])
    return similarity

def mean(s):
    total = 0.0
    for x in s:
        total += x
    return total / len(s)

def sd(s):
    m = mean(s)
    ans = 0.0
    for x in s:
        ans += (x-m) * (x-m)
    ans /= len(s)
    return ans ** 0.5


def p_value_z(randoms, actual):
    the_sd = sd(randoms)
    the_mean = mean(randoms)
    return (actual - the_mean) / the_sd


def normalize_z(scores):
    the_mean = mean(scores)
    the_sd = sd(scores)
    return [normalize_to_unit(score, the_mean, the_sd) for score in scores]

def normalize_rank(scores):
    temp = [[scores[i],i] for i in range(len(scores))]
    sorted_temp = sorted(temp, key = lambda x: x[0])
    ans = [ 0 for i in range(len(temp))]
    for i in range(len(temp)):
        ans[sorted_temp[i][1]] = float(i) / len(temp)

    return ans

def normalize_nothing(scores):
    return scores

def normalize_to_unit(score, the_mean, the_sd):
    return (score - the_mean) / the_sd


def rank(nums, x):
    count = 0
    for it in nums:
        if x > it:
            count += 1
    return float(count) / len(nums)

def get_random_weight(weights, target):
    import random
    new_weights = [0.0 for i in range(len(weights))]
    temp = range(len(weights))
    random.shuffle(temp)
    total = 0.0
    done = False
    for i in range(len(new_weights)):
        new_weights[temp[i]] = weights[temp[i]]
        total += weights[temp[i]]
        if total > target:
            break
    return new_weights


def get_alignment_mapping(seq_a, seq_b):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -0.5
    alns = pairwise2.align.globalds(seq_a, seq_b, matrix, gap_open, gap_extend)
    top_aln = alns[0]
    aln_a, aln_b, score, begin, end = top_aln
    a_i = 0
    b_i = 0
    aln_length = len(aln_a)

    a_to_b_map = {}
    b_to_a_map = {}



    for i in range(aln_length):
        if aln_a[i] != '-' and aln_b[i] != '-':
            a_to_b_map[a_i] = b_i
            b_to_a_map[b_i] = a_i
                
        if aln_a[i] != '-':
            a_i += 1
        if aln_b[i] != '-':
            b_i += 1
    return a_to_b_map, b_to_a_map


def get_random_KLs(col, seq_weights, weight_a, weight_b, pseudo_count_dict, num_trials):
    ans = []
    for i in range(num_trials):
        a = get_random_weight(seq_weights, weight_a)
        b = get_random_weight(seq_weights, weight_b)
        ans.append(get_KL_weighted(col, a, b, pseudo_count_dict, None))

    return ans

def predict_position_energy_weighted(params, mutation, use_neighbor, ignore_pos, max_neighbor, num_trials, pseudo_total, sim_f, to_neighbor_p_value):

    import wc
    import objects

    protein_name = mutation[0]
    pos = mutation[1]
    wild_res = mutation[2]
    mut_res = mutation[3]
    params.set_param('uniprot_id', protein_name)
    seq = wc.get_stuff(objects.dW, params)

    assert seq[pos] == wild_res


    import wrapper

    msa = wc.get_stuff(wrapper.my_msa_obj_wrapper, params)



    score = 0
    


    #seq_weights = [1.0 for i in range(len(msa))]

    #seq_weights = 

    #params.set_param('which_msa', 0)
    
    node_msa = wc.get_stuff(wrapper.my_msa_obj_wrapper, params)

    column = node_msa.get_column(pos)
    node_seq_weights = wc.get_stuff(objects.general_seq_weights, params)
    
    #params.set_param('which_msa', 2)
    msa = wc.get_stuff(wrapper.my_msa_obj_wrapper, params)

    neighbor_seq_weights = wc.get_stuff(objects.general_seq_weights, params)



    if not ignore_pos:
        #mut_weight = sum([seq_weights[i] for i in range(len(msa)) if msa[pos,i] == mut_res])
        #wild_weight = sum([seq_weights[i] for i in range(len(msa)) if msa[pos,i] == wild_res])

        
        #mut_count = column.count(mut_res)
        #wild_count= column.count(wild_res)
        #if wild_similarity < 1:
        #    print wild_res, mut_res, column
        #    pdb.set_trace()
        
        mut_similarity = compute_similarity_score_to_residue(column, node_seq_weights, mut_res, sim_f)
        wild_similarity = compute_similarity_score_to_residue(column, node_seq_weights, wild_res, sim_f)
        #score += math.log((mut_similarity + 1) / wild_similarity)
        score += (mut_similarity + 1) / wild_similarity

        #assert abs(score - second_score) < .001
        #score += math.log((mut_weight + 1) / (wild_weight))
    
        #score = -1.0 * mutation[-3]
    
    neighbor_score = 0
    
    if use_neighbor:

        # get neighbors/weights

        all_neighbors = wc.get_stuff(objects.general_neighbors_w_weight_w, params)
        pos_neighbors = all_neighbors[pos]
        sorted_pos_neighbors = sorted(pos_neighbors, key = lambda elt: elt[1], reverse = True)

        neighbors = [x[0] for x in sorted_pos_neighbors[0:min(max_neighbor,len(sorted_pos_neighbors))]]
        #neighbor_weights = [x[1] for x in sorted_pos_neighbors[0:min(max_neighbor,len(sorted_pos_neighbors))]]
        neighbor_weights = [1.0 for i in range(len(neighbors))]

        # get pseudo_counts
        pseudo_count_dict = {}
        for key in range(global_stuff.q):
            pseudo_count_dict[key] = pseudo_total / global_stuff.q





        # none of weights have to be normalized
        def get_neighbor_score(msa, weight_a, weight_b, neighbors, neighbor_weights, pseudo_count_dict):

            num_neighbors = len(neighbors)
            assert(len(neighbors) == len(neighbor_weights))
            score = 0
            neighbor_weights = normalize(neighbor_weights)

            #neighbors = range(num_neighbors)
            

            for i in range(num_neighbors):
                choose_neighbor_probs = {}
                for j in range(global_stuff.q):
                    choose_neighbor_probs[j] = 0
                choose_neighbor_probs[global_stuff.aa_to_num[seq[neighbors[i]]]] = 1.0
                #score += neighbor_weights[i] * get_KL_weighted(msa.get_column(neighbors[i]), weight_a, weight_b, pseudo_count_dict, choose_neighbor_probs)



                #score += neighbor_weights[i] * get_KL_weighted(msa.get_column(neighbors[i]), weight_a, weight_b, pseudo_count_dict, choose_neighbor_probs) / get_entropy_weighted(msa.get_column(neighbors[i]), neighbor_seq_weights, pseudo_count_dict)
                asdf = get_KL_weighted(msa.get_column(neighbors[i]), weight_a, weight_b, pseudo_count_dict, choose_neighbor_probs)

                random_kls = get_random_KLs(msa.get_column(neighbors[i]), neighbor_seq_weights, sum(weight_a), sum(weight_b), pseudo_count_dict, num_trials)
                #score += asdf * neighbor_weights[i]
                score += p_value_z(random_kls, asdf) * neighbor_weights[i]
                #score += neighbor_weights[i] * (asdf / mean(random_kls))
                #score += rank(random_kls, asdf) * neighbor_weights[i]
                
            return score

        
                    


        actual_weight_a = [neighbor_seq_weights[i] if msa[i,pos] == wild_res else 0.0 for i in range(len(msa))]
        actual_weight_b = [neighbor_seq_weights[i] if msa[i,pos] == mut_res else 0.0 for i in range(len(msa))]


        neighbor_cols = [ [msa[j,i] for j in range(len(msa))] for i in neighbors]


        actual_neighbor_score = get_neighbor_score(msa, actual_weight_a, actual_weight_b, neighbors, neighbor_weights, pseudo_count_dict)

        if actual_neighbor_score < 0:
            print actual_neighbor_score

        



        if to_neighbor_p_value:





            mut_weight = sum([neighbor_seq_weights[i] for i in range(len(msa)) if msa[i,pos] == mut_res])
            wild_weight = sum([neighbor_seq_weights[i] for i in range(len(msa)) if msa[i,pos] == wild_res])
            random_scores = []
            for i in range(num_trials):
                random_weight_a = get_random_weight(neighbor_seq_weights, wild_weight)
                random_weight_b = get_random_weight(neighbor_seq_weights, mut_weight)





                a_random_score = get_neighbor_score(msa, random_weight_a, random_weight_b, neighbors, neighbor_weights, pseudo_count_dict)

                random_scores.append(a_random_score)


            normalize_neighbor_by_z = True
        
            if normalize_neighbor_by_z:
                random_mean = mean(random_scores)
                random_sd = sd(random_scores)
        
                try:
                    neighbor_score = normalize_to_unit(actual_neighbor_score, random_mean, random_sd)
                except:

                    asdf=2
                    neighbor_score = 0


            else:

                neighbor_score = rank(random_scores, actual_neighbor_score)



        else:
            neighbor_score = actual_neighbor_score

    print >> sys.stderr, score, neighbor_score, len(msa)


    return (score - neighbor_score) * -1.0





def get_KL_real(d1, d2, weights):

    # keep dictionary of counts for each distribution


    pseudo_count = 1.0 / 21000


    
    d1_dict = {}
    d1_weight = 0.0
    d2_dict = {}
    d2_weight = 0.0

    for key in global_stuff.aa_to_num.keys():
        d1_dict[key] = pseudo_count
        d1_weight += pseudo_count
        d2_dict[key] = pseudo_count
        d2_weight += pseudo_count
    
    for i in range(len(d1)):
        if d1[i] not in global_stuff.ignore_aas:
            d1_weight += weights[i]
            if d1[i] in d1_dict.keys():
                d1_dict[d1[i]] = d1_dict[d1[i]] + weights[i]
            else:
                d1_dict[d1[i]] = weights[i]
                
    for k in d1_dict.keys():
        d1_dict[k] = float(d1_dict[k]) / d1_weight
        

    for i in range(len(d2)):
        if d2[i] not in global_stuff.ignore_aas:
            d2_weight += weights[i]
            if d2[i] in d2_dict.keys():
                d2_dict[d2[i]] = d2_dict[d2[i]] + weights[i]
            else:
                d2_dict[d2[i]] = weights[i]



    for k in d2_dict.keys():
        d2_dict[k] = d2_dict[k] / d2_weight

    ans = 0
    for k in global_stuff.aa_to_num.keys():
        ans += d1_dict[k] * math.log(d1_dict[k] / d2_dict[k])

    return ans
                                                                                                                                                                                                                                                        
            
def filter_column_pair(x, y, escape):

    x_skip = ''
    y_skip = ''
    assert len(x) == len(y)
    for i in range(len(x)):
        if x[i] not in escape and y[i] not in escape:
            x_skip += x[i]
            y_skip += y[i]
    return x_skip, y_skip

def read_hp_values(file_location):

    # read folder_name
    f = open(constants.INFO_FOLDER + file_location)

    hp_values = param.param()
    folder_name = f.readline().strip()



    for line in f.readlines():
        s = line.strip().split(',')
        if s[0] == 'hp':
            to_add = []
            the_type = s[2]
            the_name = s[1]
            for i in range(3,len(s)):
                if the_type == 'f':
                    to_add.append(float(s[i]))
                elif the_type == 'i':
                    to_add.append(int(s[i]))
                elif the_type == 's':
                    to_add.append(s[i])
            hp_values.set_param(the_name, to_add)

    return hp_values

def get_aux_folder(pdb_name, chain_letter, start, end):
    return constants.AUX_FOLDER + string.lower(pdb_name) + '_' + string.upper(chain_letter) + '_' + str(start) + '_' + str(end) + '/'

# returns mat normalized by columns
def normalize_mat(mat):
    height = len(mat)
    width = len(mat[0])
    normalized = [ [0 for i in range(width)] for j in range(height)]

    for i in range(width):
        a_mean_sum = 0.0
        for j in range(height):
            a_mean_sum = a_mean_sum + mat[j][i]
        a_mean = a_mean_sum / height
        a_var_sum = 0
        for j in range(height):
            a_var_sum = a_var_sum + pow((mat[j][i] - a_mean), 2)
        sd = math.sqrt(a_var_sum / height)
        if abs(sd) > .00001:
            for j in range(height):
                normalized[j][i] = (mat[j][i] - a_mean) / sd
    return normalized


def shorten(x):
    x = re.sub(r'\'','',x)
    x = re.sub(r'class','',x)
    x = re.sub(r' ','',x)
    x = re.sub(r'<','',x)
    x = re.sub(r'>','',x)
    x = re.sub(r'f\.','',x)
    x = re.sub(r'\),\(',')(',x)
    return x


def super_shorten(x):

    x = re.sub(r'\'','',x)
    x = re.sub(r'class','',x)
    x = re.sub(r' ','',x)

    
    x = re.sub(r'<','',x)
    x = re.sub(r'>','',x)
    x = re.sub(r'f\.','',x)
    x = re.sub(r'\),\(',')(',x)
    #x = re.sub(r'\)\(','|',x)
    #x = re.sub(r'\[\(','[',x)
    #x = re.sub(r'\)\]',']',x)
    #pdb.set_trace()
    return x


def get_KL(d1, d2):

    d1d = {}
    d2d = {}
    d1d2d = {}
    assert len(d1) == len(d2)
    #print >> sys.stderr, len(d1)
    for i in range(len(d1)):
        if d1[i] in d1d:
            d1d[d1[i]] += 1.0 / len(d1)
        else:
            d1d[d1[i]] = 1.0 / len(d1)
        if d2[i] in d2d:
            d2d[d2[i]] += 1.0 / len(d1)
        else:
            d2d[d2[i]] = 1.0 / len(d1)
        if (d1[i],d2[i]) in d1d2d:
            d1d2d[(d1[i],d2[i])] += 1.0 / len(d1)
        else:
            d1d2d[(d1[i],d2[i])] = 1.0 / len(d1)
    ans = 0
    for key in d1d2d:
        try:
            ans += d1d2d[key] * ((math.log(d1d2d[key]) - math.log(d1d[key[0]])) - math.log(d2d[key[1]]))
        except:
            #pdb.set_trace()
            pass
    return ans

totalcount = 0
Xcount = 0

def do_map(aa,mapping):
    try:


        return mapping[aa]
    except:

        return mapping['-']


def get_KL_fast_alt(msa, a, b, char_to_num, weights):

    

    max_num = len(char_to_num) + 1

    d1d = [0 for i in range(max_num)]
    d2d = [0 for i in range(max_num)]
    d1d2d = [ [0 for i in range(max_num)] for j in range(max_num) ]

    total = 0

    for i in range(len(msa)):
        if msa[i][a] != '_' and msa[i][b] != 'b':
            d1d[do_map(msa[i][a], global_stuff.aa_to_num)] += weights[i]
            d2d[do_map(msa[i][b], global_stuff.aa_to_num)] += weights[i]
            d1d2d[do_map(msa[i][a], global_stuff.aa_to_num)][do_map(msa[i][b], global_stuff.aa_to_num)] += weights[i]
            total += weights[i]

    for i in range(max_num):
        d1d[i] /= total
        d2d[i] /= total
        for j in range(max_num):
            d1d2d[i][j] /= total
    ans = 0

    for i in range(max_num):
        for j in range(max_num):
            if d1d2d[i][j] > 0:
                ans += d1d2d[i][j] * ( math.log(d1d2d[i][j]) - math.log(d1d[i]) - math.log(d2d[j]))
                
    return ans

def get_entropy_weighted(col, weight, pseudo_dict):

    total = 0.0
    d = {}
    for key in pseudo_dict:
        d[key] = pseudo_dict[key]
        total += pseudo_dict[key]

    for char, a_weight in zip(col, weight):
        if char not in global_stuff.ignore_aas:
            d[global_stuff.aa_to_num[char]] += a_weight
            total += a_weight

    for key in d:
        d[key] /= total

    import math
    ans = 0.0
    for key in d:
        ans += d[key] * math.log(d[key])
    return -1.0 * ans


def get_KL_weighted(col, weight_1, weight_2, pseudo_count_dict, choose_neighbor_probs):
    import pdb

    # keep dictionary of counts for each distribution
    d1_dict = {}
    d1_weight = 0.0
    d2_dict = {}
    d2_weight = 0.0
    import pdb



    def print_non_zero_weight_col(col, weight):
        assert len(col) == len(weight)
        ans = ''
        for i in range(len(col)):
            if weight[i] > 1e-10:
                ans = ans + col[i]
        print ans

    #print 'one'
    #print_non_zero_weight_col(col, weight_1)
    #print 'two'
    #print_non_zero_weight_col(col, weight_2)

    for key in pseudo_count_dict:
        d1_dict[key] = pseudo_count_dict[key]
        d1_weight += pseudo_count_dict[key]
        d2_dict[key] = pseudo_count_dict[key]
        d2_weight += pseudo_count_dict[key]
    

    

    try:
        for i in range(len(col)):

            if col[i] not in global_stuff.ignore_aas:
                d1_weight += weight_1[i]
                if global_stuff.aa_to_num[col[i]] in d1_dict.keys():
                    d1_dict[global_stuff.aa_to_num[col[i]]] = d1_dict[global_stuff.aa_to_num[col[i]]] + weight_1[i]
                else:
                    d1_dict[global_stuff.aa_to_num[col[i]]] = weight_1[i]

                d2_weight += weight_2[i]
                if global_stuff.aa_to_num[col[i]] in d2_dict.keys():
                    d2_dict[global_stuff.aa_to_num[col[i]]] = d2_dict[global_stuff.aa_to_num[col[i]]] + weight_2[i]
                else:
                    d2_dict[global_stuff.aa_to_num[col[i]]] = weight_2[i]

    except:
        raise Exception
        pdb.set_trace()




    try:
        for k in d1_dict.keys():
            d1_dict[k] = float(d1_dict[k]) / d1_weight
    except Exception, err:
        raise Exception

    
    for k in d1_dict.keys():
        d2_dict[k] = float(d2_dict[k]) / d2_weight

    


    ans = 0


    #choose_neighbor_probs = d1_dict
    
    try:
        for k in d1_dict.keys():
            #if choose_neighbor_probs[k] > 0:
            #    ans += choose_neighbor_probs[k] * math.log(d1_dict[k] / d2_dict[k])
                
            
            if d1_dict[k] > 0:
                ans += d1_dict[k] * math.log(d1_dict[k] / d2_dict[k])
    except Exception, err:
        print err
        print >> sys.stderr, d1_dict[k], d2_dict[k], k

        x=2

    if ans < 0:
        print ans
        pdb.set_trace()
        
    return ans 




def physical_distance(x):
    norm = 0;
    for i in range(len(x)):
        norm = norm + x[i] * x[i]
    return math.sqrt(norm)

def print_stuff_dec(f):

    def g(*args, **kwargs):
        #print >> sys.stderr, 'calling ', f.func_name, ' with ', args, kwargs
        ans = f(*args, **kwargs)
        #print >> sys.stderr, f.func_name, ' returned ', ans
        return ans
    
    return g

def get_object(p_wrapper, params, recalculate = False, to_pickle = True, use_pickle = True):
    return the_obj_manager.get_variable(p_wrapper(params), recalculate, to_pickle, use_pickle)

def get_file(p_wrapper, params, recalculate = False, option = 'r'):
    return the_file_manager.get_file(p_wrapper(params), recalculate, option)

def get_transpose(mat):
    height = len(mat)
    width = len(mat[0])
    m = [ [-1 for i in range(height)] for j in range(width)]
    for i in range(width):
        for j in range(height):
            m[i][j] = mat[j][i]
    return m

def dict_deep_copy(d):
    to_return = {}
    for key in d.keys():
        to_return[key] = d[key]
    return to_return

def list_union(a, b):
    A = set(a)
    B = set(b)
    return list(A - B)


def write_mat(mat, f_name, the_sep = ',', option = 'w'):
    f = open(f_name, option)
    #print >> sys.stderr, mat

    for row in mat:
        
        line = string.join([('%.5f' % x)  for x in row], sep=the_sep)
        line = line + '\n'
        f.write(line)
    f.close()

def write_mat_raw(mat, f_name, the_sep = ',', option = 'w'):
    f = open(f_name, option)
    #print >> sys.stderr, mat
    for row in mat:
        
        line = string.join([str(x)  for x in row], sep=the_sep)
        line = line + '\n'
        f.write(line)
    f.close()


def read_mat_to_int_float_tuple(f):
    import pdb

    f = open(f.name, 'r')
    ans = []
    for line in f:
        this = []
        if len(line.strip()) > 0:
            try:
                s = line.strip().split(',')
                for it in s:
                    sp = it[1:len(it)-1]
                    spp = sp.split('-')
                    a = int(spp[0])
                    b = float(spp[1])
                    this.append((a,b))
            except Exception, err:
                print >> sys.stderr, err

                print >> sys.stderr, s
                
        ans.append(this)
    return ans

def write_int_float_tuple_mat(mat, f_name):
    def g(t):
        return '(' + str(t[0]) + '-' + ('%.3f' % t[1]) + ')'
    f = open(f_name, 'w')
    for row in mat:
        line = string.join([g(x) for x in row],sep = ',')
        line = line + '\n'
        f.write(line)
    f.close()
        
def write_iterable_vert(obj, f):
    f = open(f,'w')
    for x in obj:
        f.write(str(x) + '\n')
    f.close()

def write_vect(vect, f_name, the_sep = ',', option = 'w'):
    f = open(f_name, option)
    line = string.join([str(x) for x in vect], sep=the_sep)
    f.write(line)
    f.close()

def read_edge_to_int(f):
    f = open(f.name)
    ans = {}
    for line in f:
        s = line.strip().split(',')
        u = int(s[0])
        v = int(s[1])
        val = int(s[2])
        ans[(u,v)] = val
    return ans

def write_edge_to_int(obj, f_name):
    f = open(f_name, 'w')
    for key in obj:
        f.write(str(key[0]) + ',' + str(key[1]) + ',' + str(obj[key]) + '\n')
    f.close()

def read_vect_to_float(f, the_sep = ','):
    r = csv.reader(f, delimiter = the_sep)
    line = r.next()
    vect = [float(x) for x in line]
    f.close()
    return vect

def read_mat_to_float(f, the_sep = ','):
    r = csv.reader(f, delimiter = the_sep)
    mat = []
    for line in r:
        vect = [float(x) for x in line]
        mat.append(vect)
    f.close()
    return mat

def read_vect_to_int(f, the_sep = ','):
    r = csv.reader(f, delimiter = the_sep)
    line = r.next()
    vect = [int(x) for x in line]
    f.close()
    return vect

def read_mat_to_int(f, the_sep = ','):
    r = csv.reader(f, delimiter = the_sep)
    mat = []
    for line in r:
        vect = [int(x) for x in line]
        mat.append(vect)
    f.close()
    return mat

def get_overlap(n1, n2):
    n1set = set()
    n2set = set()
    count = 0
    for i in range(len(n1)):
        for it in n1[i]:
            n1set.add((i,it[0]))
            count += 1
    for i in range(len(n2)):
        for it in n2[i]:
            n2set.add((i,it[0]))
    intersect = n1set & n2set

    return len(intersect), count

def get_file_string_set(f_name):
    f = open(f_name, 'r')
    ans = set()
    for line in f:
        ans.add(line.strip())
    f.close()
    return ans

def get_representative_atom(res):
    if 'CA' in res.child_dict.keys():
        return res['CA']
    elif 'CB' in res.child_dict.keys():
        return res['CB']
    else:
        print >> sys.stderr, 'no CA or CB atom in residue'
            #pdb.set_trace()
        return res.child_list[0]


class file_sender(object):

    def __init__(self, lock_file, buildup_size):
        self.lock_file = lock_file
        self.buildup_size = buildup_size
        self.buildup = []

    def __send(self, it):
        # first try to make the remote directory
        
        import ssh
        
        here_file = it[0]
        there_file = it[1]
        there_folder = it[2]
        hostname = it[3]
        username = it[4]
        password = it[5]
        port = it[6]
        wrapper = it[7]
        object_key = it[8]
        to_remove = it[9]

        import pdb
        
        
        """
        client = ssh.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(ssh.AutoAddPolicy())
        client.connect(hostname, port, username)
        client.exec_command('mkdir ' + there_folder)
        
        cmd = 'scp ' + '\''+here_file+'\'' + ' ' + username + '@' + hostname + ':' + '\''+there_file+'\''
        import subprocess
        subprocess.call(cmd,shell=True,executable='/bin/bash')

        """
        key = ssh.RSAKey(filename='/home/fw27/.ssh/id_rsa')
        t = ssh.Transport((hostname, port))
        t.connect(username=username)
        t.auth_publickey('fultonw',key)

        sftp = ssh.SFTPClient.from_transport(t)
        try:
            print >> sys.stderr, '\t\t\tsending:', here_file, there_file
            sftp.put(here_file, there_file)
            wrapper.record_transferred_file(object_key)
        except Exception, err:
            print >> sys.stderr, err
            print >> sys.stderr, '\t\t\tfailed to send:', here_file, there_file

        if to_remove:
            try:
                import subprocess
                print >> sys.stderr, '               removing:', here_file
                subprocess.call(['rm', here_file])

            except Exception, err:
                print >> sys.stderr, err
        

        try:
            sftp.close()
        except Exception, err:
            print err
        try:
            t.close()
        except Exception, err:
            print err
        
        
        
        
    def send(self, here_file, there_file, hostname, there_folder, username, password, port, wrapper, object_key, whether_to_delete):
        #first check if if it is a file
        import pdb

        import os
        if os.path.isfile(here_file):
            import pdb

            self.buildup.append([here_file, there_file, there_folder, hostname, username, password, port, wrapper, object_key, wrapper, object_key, whether_to_delete])

        import wc
        dist_count = 0
        for it in self.buildup:
            import objects
            if type(it[0]) == objects.general_distance:
                dist_count += 1
                
        if dist_count > 0 or len(self.buildup) > self.buildup_size:
        #if len(self.buildup) > self.buildup_size:
            import FileLock
            print >> sys.stderr, "trying to send before lock: ", here_file
            with FileLock.FileLock(self.lock_file, timeout=200) as lock:
                for it in self.buildup:
                    self.__send(it)
            self.buildup = []
            import pdb

def hamming_distance(s1, s2):
    assert(len(s1) == len(s2))
    count = 0.0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            count += 1
    return count / len(s1)


def get_num_wild(msa):
    ans = [ 0 for i in range(msa.get_alignment_length())]
    for i in range(msa.get_alignment_length()):
        c = msa.get_column(i)
        ans[i] = c.count('-')
    return ans



# works with unimputed msa
def filter_msa(msa, cut_off):

    # choose a random order to add sequences
    order = random.shuffle(range(len(msa)))
    sequences = []
    added_indicies = []
    # convert msa to list of sequences
    new_msa = [ [msa[i][j] for j in range(msa.get_alignment_length())] for i in range(len(msa))]

    for i in range(len(msa)):
        to_add = True
        candidate = new_msa[i]

        if msa[i].id == 'QUERY':
            to_add = True
        else:
            for j in range(len(added_indicies)):
                if hamming_distance(candidate, new_msa[added_indicies[j]]) < cut_off:
                    to_add = False
                    break
        if to_add:
            sequences.append(candidate)
            added_indicies.append(i)
    sequences = []
    for idx in added_indicies:
        sequences.append(msa[idx])

    return MultipleSeqAlignment(sequences)


# for cosmic, mutation is of form name, pos, wild, mut, # of that mutation, # of that site, # of that gene

def mutation_to_class(mutation):
    asdf = int( mutation[8] == 1 or mutation[9] == 1 or mutation[10] == 1)

    return [asdf,int(mutation[4]>1), mutation[8], mutation[9], mutation[10]]

def cbs_to_class(mutation):
    if mutation[4] < 3.0:
        return [1]
    else:
        return [0]

def temp_cosmic(mutation):
    asdf = int( mutation[8] == 1 or mutation[9] == 1 or mutation[10] == 1)

    return [asdf,int(mutation[4]>1), mutation[8], mutation[9], mutation[10]]
    return [mutation[4], mutation[8], mutation[9], mutation[10]]
    if mutation[4] > 1:
        return 1
    else:
        return 0

def p53_stone_to_class(mutation):
    if mutation[4] >= 5:
        return [1]
    else:
        return [0]

def lys_stone_to_class(mutation):
    if mutation[4] == 2:# or mutation[4] == 1:
        return [0]
    else:
        return [1]

def protease_stone_to_class(mutation):
    if mutation[4] == 'Positive' or mutation[4] == 'Intermediate':
        return [1]
    else:
        return [0]

def reverse_stone_to_class(mutation):
    if mutation[4] == 0 or mutation[4] == 1:
        return [0]
    else:
        return [1]

def hemo_stone_to_class(mutation):
    if mutation[4] == 'Normal':
        return [0]
    else:
        return [1]

def saapdb_to_class(mutation):
    return mutation[4]

def p53_to_class(mutation):
    change = sum(mutation[4:]) / 8.0
    print 'CHANGE: ', change
    if change <= 20.0:
        return [1]
    else:
        return [0]

import numpy

class my_msa:

    # mat is list of strings where each string is a column
    def __init__(self, mat):
        self.mat = mat

    @classmethod
    def msa_from_file(cls, file):

        f = open(file.name, 'r')
        first = f.readline()
        s = first.split(' ')
        length = int(s[0])
        alignment_length = int(s[1])
        mat = [None for j in range(alignment_length)]
        i = 0
        for line in f:
            mat[i] = line.strip()
            assert len(mat[i]) == length
            i += 1
        
        return cls(mat)

    def __len__(self):
        try:
            return len(self.mat[0])
        except:
            return 0

    def __getitem__(self, ij):

        return self.mat[ij[1]][ij[0]]


    def get_alignment_length(self):
        return len(self.mat)

    def write_to_file(self, file):
        f = open(file, 'w')
        f.write(str(self.__len__()) + ' ' + str(self.get_alignment_length()) + '\n')
        for column in self.mat:
            f.write(column + '\n')
        f.close()

    def get_column(self, i):
        return self.mat[i]
                    
def parse_p_input(p, arg_string):

    for z in range(0,len(arg_string),3):
        param_name = arg_string[z]
        param_type = arg_string[z+2]
        val = arg_string[z+1]
        if param_type == 'i':
            val = int(val)
        elif param_type == 'f':
            val = float(val)
        elif param_type == 's':
            val = val
        p.set_param(param_name, val)

def conv_seq(in_file, out_file, format):
    subprocess.call(global_stuff.CONVSEQ_PATH + ' ' + in_file + ' ' + out_file + ' ' + format, shell=True, executable='/bin/bash')

def cp_file(from_file, to_file):
    subprocess.call('cp ' + '\''+from_file+'\' ' + to_file, shell=True, executable='/bin/bash')

def rm_files(the_files):
    for the_file in the_files:
        subprocess.call('rm ' + '\''+the_file+'\'', shell=True, executable='/bin/bash')
    
