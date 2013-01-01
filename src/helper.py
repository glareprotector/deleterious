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

from Bio.Align import MultipleSeqAlignment

# also sets global_stuff.RESULTS_FOLDER to proper value
def read_param(file_location):
    
    # read folder_name
    f = open(constants.INFO_FOLDER + file_location)
    the_params = param.param({})
    hp_values = param.param()
    folder_name = f.readline().strip()

    global_stuff.RESULTS_FOLDER = global_stuff.RESULTS_BASE_FOLDER + folder_name + '/'

    for line in f.readlines():
        print line
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
        pdb.set_trace()
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
        print score
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


def normalize_to_unit(score, the_mean, the_sd):
    return (score - the_mean) / the_sd


def rank(nums, x):
    count = 0
    for it in nums:
        if x > it:
            count += 1
    return float(count) / len(nums)


def predict_position_energy_weighted(params, recalculate, mutation, use_neighbor, ignore_pos, max_neighbor, weighted, num_trials, pseudo_total, sim_f):

    import wc
    import objects

    protein_name = mutation[0]
    pos = mutation[1]
    wild_res = mutation[2]
    mut_res = mutation[3]
    params.set_param('uniprot_id', protein_name)
    seq = wc.get_stuff(objects.dW, params, recalculate, False, False)

    msa = wc.get_stuff(objects.general_msa, params, recalculate, False, False)

    score = 0
    
    seq_weights = wc.get_stuff(objects.general_seq_weights, params, False, False, False)
    #seq_weights = [1.0/len(msa) for i in range(len(msa))]

    column = msa.get_column(pos)
    
    if not ignore_pos:
        #mut_weight = sum([seq_weights[i] for i in range(len(msa)) if msa[i][pos] == mut_res])
        #wild_weight = sum([seq_weights[i] for i in range(len(msa)) if msa[i][pos] == wild_res])
        mut_similarity = compute_similarity_score_to_residue(column, seq_weights, mut_res, sim_f)
        wild_similarity = compute_similarity_score_to_residue(column, seq_weights, wild_res, sim_f)
        score += math.log((mut_similarity + 1) / wild_similarity)
        #score += math.log((mut_weight + 1) / (wild_weight))
    
    return score
    #return 0.0 - score
    actual_neighbor_score = 0
    
    if use_neighbor:

        # get neighbors/weights
        all_neighbors = wc.get_stuff(objects.neighbors_w_weight_w, params, recalculate, False, False)
        pos_neighbors = all_neighbors[pos]
        sorted_pos_neighbors = sorted(pos_neighbors, key = lambda elt: elt[1], reverse = True)
        neighbors = [x[0] for x in sorted_pos_neighbors[0:min(max_neighbor,len(sorted_pos_neighbors))]]
        neighbor_weights = [x[1] for x in sorted_pos_neighbors[0:min(max_neighbor,len(sorted_pos_neighbors))]]

        # get pseudo_counts
        pseudo_count_dict = {}
        for key in global_stuff.aa_to_num:
            pseudo_count_dict[key] = pseudo_total / global_stuff.q

        # none of weights have to be normalized
        def get_neighbor_score(msa, weight_a, weight_b, neighbors, neighbor_weights, pseudo_count_dict):
            num_neighbors = len(neighbors)
            assert(len(neighbors) == len(neighbor_weights))
            score = 0
            neighbor_weights = normalize(neighbor_weights)

            for i in range(num_neighbors):
                score += neighbor_weights[i] * get_KL_weighted(msa, weight_a, weight_b, neighbors[i], pseudo_count_dict)
            return score

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
                    


        actual_weight_a = [seq_weights[i] if msa[i][pos] == wild_res else 0.0 for i in range(len(msa))]
        actual_weight_b = [seq_weights[i] if msa[i][pos] == mut_res else 0.0 for i in range(len(msa))]

        actual_neighbor_score = get_neighbor_score(msa, actual_weight_a, actual_weight_b, neighbors, neighbor_weights, pseudo_count_dict)
        #return 0 - actual_neighbor_score
        #return score - actual_neighbor_score

        """
        constraints_a = [(pos,wild_res)]
        filter_a_msa = filter_msa_based_on_pos_constraint(msa, constraints_a)
        constraints_b = [(pos,mut_res)]
        filter_b_msa = filter_msa_based_on_pos_constraint(msa, constraints_b)
        #all_neighbors = wc.get_stuff(objects.neighbors_w, params, recalculate, False, False)
        neighbor = neighbors[0]
        na_col = filter_a_msa.get_column(neighbor)
        nb_col = filter_b_msa.get_column(neighbor)
        na_col_no_skip = [x for x in na_col if x != '-']
        nb_col_no_skip = [x for x in nb_col if x != '-']
        asdf = get_KL_real(nb_col_no_skip, na_col_no_skip, seq_weights)
        """


        



        #return actual_score
        mut_weight = sum([seq_weights[i] for i in range(len(msa)) if msa[i][pos] == mut_res])
        wild_weight = sum([seq_weights[i] for i in range(len(msa)) if msa[i][pos] == wild_res])
        random_scores = []
        for i in range(num_trials):
            random_weight_a = get_random_weight(seq_weights, wild_weight)
            random_weight_b = get_random_weight(seq_weights, mut_weight)

            a_random_score = get_neighbor_score(msa, random_weight_a, random_weight_b, neighbors, neighbor_weights, pseudo_count_dict)

            random_scores.append(a_random_score)
        """
        ans = rank(random_scores, actual_neighbor_score)
        print ans
        return ans
        """


        

        random_mean = mean(random_scores)
        random_sd = sd(random_scores)
        
        try:
            neighbor_score = normalize_to_unit(actual_neighbor_score, random_mean, random_sd)
        except:
            pdb.set_trace()
            asdf=2
            neighbor_score = 0


        normalized_randoms = [normalize_to_unit(x, random_mean, random_sd) for x in random_scores]



#        x = 0
#        for score in random_scores:
#            if actual_score > score:
#                x += 1
#        neighbor_score = float(x) / num_trials
        print neighbor_score, len(msa)

        return score - neighbor_score


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
    #print len(d1)
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

def get_KL_weighted(msa, weight_1, weight_2, neighbor, pseudo_count_dict):

    # keep dictionary of counts for each distribution
    d1_dict = {}
    d1_weight = 0.0
    d2_dict = {}
    d2_weight = 0.0
    import pdb

    """
    for key in pseudo_count_dict:
        d1_dict[key] = pseudo_count_dict[key]
        d1_weight += pseudo_count_dict[key]
        d2_dict[key] = pseudo_count_dict[key]
        d2_weight += pseudo_count_dict[key]
    """

    

        
    for i in range(len(msa)):
        if msa[i][neighbor] not in global_stuff.ignore_aas:
            d1_weight += weight_1[i]
            if global_stuff.aa_to_num[msa[i][neighbor]] in d1_dict.keys():
                d1_dict[global_stuff.aa_to_num[msa[i][neighbor]]] = d1_dict[global_stuff.aa_to_num[msa[i][neighbor]]] + weight_1[i]
            else:
                d1_dict[global_stuff.aa_to_num[msa[i][neighbor]]] = weight_1[i]

            d2_weight += weight_2[i]
            if global_stuff.aa_to_num[msa[i][neighbor]] in d2_dict.keys():
                d2_dict[global_stuff.aa_to_num[msa[i][neighbor]]] = d2_dict[global_stuff.aa_to_num[msa[i][neighbor]]] + weight_2[i]
            else:
                d2_dict[global_stuff.aa_to_num[msa[i][neighbor]]] = weight_2[i]

#    pseudo_count = 1.0 / len(msa)
    #pseudo_count = 1.0
#    all_keys = set(d1_dict.keys()) | set(d2_dict.keys())
#    for key in d1_dict.keys():
#        d2_dict[key] += pseudo_count
#        d2_weight += pseudo_count


    try:
        for k in d1_dict.keys():
            d1_dict[k] = float(d1_dict[k]) / d1_weight
    except Exception, err:
        raise Exception

    
    for k in d1_dict.keys():
        d2_dict[k] = float(d2_dict[k]) / d2_weight

    


    ans = 0
    try:
        for k in d1_dict.keys():
            if d1_dict[k] > 0:
                ans += d1_dict[k] * math.log(d1_dict[k] / d2_dict[k])
    except:
        print d1_dict[k], d2_dict[k], k
        pdb.set_trace()
        x=2
    return ans




def physical_distance(x):
    norm = 0;
    for i in range(len(x)):
        norm = norm + x[i] * x[i]
    return math.sqrt(norm)

def print_stuff_dec(f):

    def g(*args, **kwargs):
        #print 'calling ', f.func_name, ' with ', args, kwargs
        ans = f(*args, **kwargs)
        #print f.func_name, ' returned ', ans
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
    #print mat
    for row in mat:
        
        line = string.join([('%.5f' % x)  for x in row], sep=the_sep)
        line = line + '\n'
        f.write(line)
    f.close()

def write_mat_raw(mat, f_name, the_sep = ',', option = 'w'):
    f = open(f_name, option)
    #print mat
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
                print err
                pdb.set_trace()
                print s
                
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
        print 'no CA or CB atom in residue'
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



        client = ssh.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(ssh.AutoAddPolicy())
        client.connect(hostname, port, username, password)
        client.exec_command('mkdir ' + there_folder)


        t = ssh.Transport((hostname, port))
        t.connect(username=username, password=password)

        sftp = ssh.SFTPClient.from_transport(t)
        try:
            print '\t\t\tsending:', here_file, there_file
            sftp.put(here_file, there_file)
            wrapper.record_transferred_file(object_key)
        except Exception, err:
            print err
            print '\t\t\tfailed to send:', here_file, there_file

        try:
            import subprocess
            print '               removing:', here_file
            subprocess.call(['rm', here_file])

        except Exception, err:
            print err


        
        
    def send(self, here_file, there_file, hostname, there_folder, username, password, port, wrapper, object_key):
        #first check if if it is a file
        import pdb

        import os
        if os.path.isfile(here_file):
            import pdb

            self.buildup.append([here_file, there_file, there_folder, hostname, username, password, port, wrapper, object_key, wrapper, object_key])

        import wc
        dist_count = 0
        for it in self.buildup:
            import objects
            if type(it[0]) == objects.general_distance:
                dist_count += 1
                
        if dist_count > 0 or len(self.buildup) > self.buildup_size:
        #if len(self.buildup) > self.buildup_size:
            import FileLock
            print "trying to send before lock: ", here_file
            with FileLock.FileLock(self.lock_file, timeout=2) as lock:
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
    return mutation[4]

def temp_cosmic(mutation):
    if mutation[5] > 1:
        return 1
    else:
        return 0



                    
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
