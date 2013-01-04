
from wrapper_decorator import dec
import pdb
import wrapper
import param
import numpy




import Bio.PDB

import global_stuff
import helper
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.PDB import Polypeptide
from Bio.Blast import NCBIXML
from Bio.Align.Applications import MuscleCommandline
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

import math
import subprocess
import string
import os
import random
import pdb
import re


import global_stuff



class bW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['uniprot_id'])

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        # all this will do is read specified fasta file, and save it in my format
        protein_name = self.get_param(params, 'uniprot_id')
        protein_folder = global_stuff.base_folder + protein_name + '/'
        existing_seq_file = protein_folder + 'seq'

        old_seq = SeqIO.read(open(existing_seq_file,'r'), 'fasta')

        seq_record = Bio.SeqRecord.SeqRecord(old_seq)
        SeqIO.write(old_seq, self.get_holding_location(), 'fasta')
        return open(self.get_holding_location(),'r')

class dW(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return bW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        seq_file_handle = self.get_var_or_file(bW, params, recalculate, True, False)
        asdf = SeqIO.read(seq_file_handle, 'fasta')
        if asdf[-1] == '*':
            asdf = asdf[0:len(asdf)-1]
        return asdf

    def whether_to_override(self, object_key):
        return False

# blast results file wrapper(xml format)
class adW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['ev', 'which_blast']) | bW.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return False
        #if the file size is too small, we know there was something wrong
        import os
        location = self.get_file_location(object_key)
        if os.path.getsize(location) < 1:
            return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        seq_records = []
        f = self.get_var_or_file(bW, params, recalculate, False, False, False)
        query = SeqIO.parse(f, 'fasta')
        seq_records.append(query)

        print 'RUNNING BLAST!!!!!!!'
        if self.get_param(params, 'which_blast') == 0:
            psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLAST_PATH, outfmt = 5, query = '\''+f.name+'\'', db = global_stuff.BLASTDB_PATH, out = self.get_holding_location(), evalue = self.get_param(params, 'ev'), num_iterations = self.get_param(params, 'psiblast_iter'))
        elif self.get_param(params, 'which_blast') == 1:
            #psi_blast_cline = global_stuff.BLASTP_PATH + 
            psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLASTP_PATH, outfmt = 5, query = '\''+f.name+'\'', db = global_stuff.BLASTDB_PATH, out = self.get_holding_location(), evalue = self.get_param(params, 'ev'))
        #psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLAST_PATH, outfmt = 5, query = '\''+f.name+'\'', out = self.get_holding_location())
        #pdb.set_trace()
        print psi_blast_cline
        subprocess.call(str(psi_blast_cline), shell=True, executable='/bin/bash')

        return open(self.get_holding_location())


# processsed blast results format(bunch of sequences in fasta in 1 file)
class aeW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        keys = adW.get_all_keys(params, self) | dW.get_all_keys(params, self) | set(['blmax', 'ev']) 
        return keys

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        # parse blast xml file, then do processing
        
        blast_xml_handle = self.get_var_or_file(adW, params, recalculate, False, False, False)
        if self.get_param(params, 'which_blast') == 1:
            record = NCBIXML.read(blast_xml_handle)
        elif self.get_param(params, 'which_blast') == 0:
            records = NCBIXML.parse(blast_xml_handle)
            iters = [x for x in records]
            record = iters[-1]
        seen = set()
        seq_records = []
        # add the query sequence, and have a set so that only add each sequence once

        query = self.get_var_or_file(dW, params, recalculate, True, False, False)
        query.id = 'QUERY'
        seq_records.append(query)
        seen.add(query.seq.tostring())
        # add high scoring pairs in alignments with sufficiently low evalue that have not been seen
        # add them to temp data structure so that i can sort by evalue before adding them
        temp = []
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < self.get_param(params, 'ev') and not hsp.sbjct in seen:
                    temp.append([hsp.expect, Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(hsp.sbjct), id = alignment.hit_id)])
        temp_sorted = sorted(temp, key = lambda x:x[0])
        for i in range(min(self.get_param(params, 'blmax'), len(temp_sorted))):
            seq_records.append(temp_sorted[i][1])
                
        # write hits to fasta file
        output_handle = open(self.get_holding_location(), 'w')
        SeqIO.write(seq_records, output_handle, 'fasta')
        print 'WROTE ', self.get_holding_location()

        return output_handle


# gets the result of muscle
class afW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return aeW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        msa_input_handle = self.get_var_or_file(aeW, params, recalculate, to_pickle, to_filelize, always_recalculate)
        cline = MuscleCommandline(cmd = global_stuff.MUSCLE_PATH, input = '\''+msa_input_handle.name+'\'', out = self.get_holding_location(), clw = False, maxiters = 2)
        #cline()

        subprocess.call(str(cline), shell=True, executable='/bin/bash')
        return open(self.get_holding_location())




class psicov_input_file(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return general_msa.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        msa = self.get_var_or_file(general_msa, params, False, False, False)
        # first, find the query index
        query_idx = None
        for i in range(len(msa)):
            if msa[i].name == 'QUERY':
                query_idx = i
                break
        assert query_idx != None
        f = open(self.get_holding_location(), 'w')
        f.write(msa[query_idx].seq.tostring() + '\n')
        for i in range(len(msa)):
            if i != query_idx:
                f.write(msa[i].seq.tostring() + '\n')
        f.close()
        return open(self.get_holding_location(), 'r')

    def whether_to_override(self, object_key):
        return True

class psicov_output_file(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return psicov_input_file.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        separating_dist = self.get_param(params, 'psicov_sep')
        gap_ignore = self.get_param(params, 'psicov_gap')
        r = self.get_param(params, 'psicov_r')
        input_file = self.get_var_or_file(psicov_input_file, params, False, False, False)
        cmd = global_stuff.PSICOV_PATH + ' -p ' + ' -r ' + str(r) + ' -j ' + str(separating_dist) + ' -g ' + str(gap_ignore) + ' ' + '\''+input_file.name+'\'' + ' > ' + self.get_holding_location()
        pdb.set_trace()
        subprocess.call(cmd, shell=True, executable='/bin/bash')
        return open(self.get_holding_location(), 'r')

    def whether_to_override(self, object_key):
        return True

class psicov_distance(wrapper.mat_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return psicov_output_file.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        psicov_output = self.get_var_or_file(psicov_output_file, params, False, False, False)
        f = open(psicov_output.name, 'r')
        seq = self.get_var_or_file(dW, params, False, False, False)
        dist = [ [0.0 for i in range(len(seq))] for j in range(len(seq))]
        for line in f:
            s = line.strip().split(' ')
            i = int(s[0])
            j = int(s[1])
            val = float(s[4])
            dist[i][j] = val
            dist[j][i] = val
        return dist


    def whether_to_override(self, object_key):
        return True

class msa_file(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return general_msa.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        msa = self.get_var_or_file(general_msa, params, False, False, False)
        f = open(self.get_holding_location(),'w')
        for i in range(len(msa)):
            f.write(msa[i].seq.tostring()+'\n')
        f.close()
        #AlignIO.write(msa, self.get_holding_location(), 'fasta')
        return open(self.get_holding_location(),'r')

class MIP_input_msa(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['which_impute', 'filter_co']) | general_msa.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        # will read in each sequence in msa, and write to new msa, except that if the name is query, will write it with a specified GI number instead
        # will also replace gaps
        msa = self.get_var_or_file(general_msa, params, False, False, False)

        filtered_msa = helper.filter_msa(msa, self.get_param(params, 'filter_co'))

        sequences = []

        which_impute = self.get_param(params, 'which_impute')
        for seq in filtered_msa:
            s = seq.seq.tostring()
            new_s = ''
            if which_impute == 0:
                for i in range(len(s)):
                    if s[i] in global_stuff.ignore_aas:
                        new_s += global_stuff.MIP_wild_char
                    else:
                        new_s += s[i]
            elif which_impute == 1:
                for i in range(len(s)):
                    if s[i] in global_stuff.ignore_aas:
                        new_s += random.choice(global_stuff.aa_to_num.keys())
                    else:
                        new_s += s[i]

            assert new_s.count('X')==0
            if seq.name != 'QUERY':
                the_id = seq.id
            else:

                the_id = 'gi|' + global_stuff.query_gi_number + '|a|a|a|a'
            sequences.append(SeqIO.SeqRecord(seq=SeqIO.Seq(new_s), id = the_id))

        return MultipleSeqAlignment(sequences)

    def whether_to_override(self, object_key):
        return True

class MIP_input_msa_file(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return MIP_input_msa.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        msa = self.get_var_or_file(MIP_input_msa, params, False, False, False)
        AlignIO.write(msa, self.get_holding_location(), 'fasta')
        return open(self.get_holding_location(), 'r')

    def whether_to_override(self, object_key):
        return True

class MIP_distance_file(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return MIP_input_msa_file.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        modified_alignment_file = self.get_var_or_file(MIP_input_msa_file, params, False, False, False)
        cmd = global_stuff.MIP_PATH + ' -i ' + '\'' + modified_alignment_file.name + '\'' + ' -o ' + self.get_holding_location() + ' -n ' + str(1) + ' -g ' + global_stuff.query_gi_number
        print cmd
        subprocess.call(cmd, shell=True, executable='/bin/bash')
        subprocess.call('mv ' + self.get_holding_location() + '_MIp.txt ' + self.get_holding_location(), shell=True, executable='/bin/bash')
        subprocess.call('rm ' + self.get_holding_location() + '_count.txt', shell=True, executable='/bin/bash')
        subprocess.call('rm ' + self.get_holding_location() + '_MIp.dot', shell=True, executable='/bin/bash')
        return open(self.get_holding_location(), 'r')

    def whether_to_override(self, object_key):
        return True


class MIP_distance(wrapper.mat_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return MIP_distance_file.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        f = self.get_var_or_file(MIP_distance_file, params)
        #skip first line
        f.next()
        # get length of sequence
        seq = self.get_var_or_file(dW, params)
        length = len(seq)
        distances = [[0 for i in range(length)] for j in range(length)]
        for line in f:
            s = line.strip().split('\t')
            x = int(s[0]) - 1
            y = int(s[1]) - 1
            d = float(s[8])
            distances[x][y] = d
        return distances


# processed msa output(columns with skips removed)
class agW(wrapper.msa_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return afW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        f = self.get_var_or_file(afW, params, recalculate, False, False, False)
        msa = AlignIO.read(f.name, 'fasta')
        # search for the query sequence
        idx = -1
        for i in range(len(msa)):
            if msa[i].id == 'QUERY':
                idx = i
                break
        # find the first non-insertion column
        i = 0

        while msa[idx,i] == '-':
            #print msa[idx,i]
            i = i + 1
            #pdb.set_trace()

        to_return = msa[:,i:(i+1)]
        # add in all the other columns

        for k in range(i+1, msa.get_alignment_length()):
            if msa[idx,k] != '-':
                #print k
                to_return = to_return + msa[:,k:(k+1)]

        return to_return


class their_agW(wrapper.msa_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['uniprot_id'])

    def whether_to_override(self, object_key):
        return False

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
      
        name = self.get_param(params, 'uniprot_id')
        their_file = global_stuff.polyphen_msa_directory + name + '.aln.mine'
        #their_file = global_stuff.base_folder + name + '/fake_msa_short'
        f = open(their_file, 'r')
        msa = AlignIO.read(f, 'fasta')
        return msa

class easy_to_read_msa_format(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return agW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        f = open(self.get_holding_location(), 'w')
        recalculate = False
        #pdb.set_trace()
        msa = self.get_var_or_file(agW, params, recalculate, True, True, False)

        f.write(str(len(msa)) + '\t' + str(msa.get_alignment_length()) + '\n')
        for i in range(len(msa)):
            to_write = ''
            for j in range(msa.get_alignment_length()):
                to_write = to_write + msa[i,j]
            #print len(to_write), to_write
            f.write(to_write + '\n')
        return f

class pairwise_dist(wrapper.mat_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return general_msa.get_all_keys(params, self) | general_seq_weights.get_all_keys(params, self)
            

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        import pdb


        msa = self.get_var_or_file(general_msa, params, False, True, False, False)
        print 'length: ', msa.get_alignment_length(), len(msa)
        import pdb
        import datetime
        start = datetime.datetime.now()
        dists = [ [0 for j in range(msa.get_alignment_length())] for i in range(msa.get_alignment_length())]

        new_msa = [ [msa[i][j] for j in range(msa.get_alignment_length())] for i in range(len(msa))]

        weights = self.get_var_or_file(general_seq_weights, params, False, False, False)



        print datetime.datetime.now() - start

        for i in range(msa.get_alignment_length()):
            for j in range(msa.get_alignment_length()):
                x_no_skip = ''
                y_no_skip = ''
                """
                for k in range(len(msa)):
                    if msa[k][i] != '-' and msa[k][j] != '-':
                        x_no_skip += msa[k][i]
                        y_no_skip += msa[k][j]
                """
                #dists[i][j] = helper.get_KL(x_no_skip, y_no_skip)
                #dists[i][j] = helper.get_KL_fast(msa.get_column(i), msa.get_column(j), global_stuff.aa_to_num)
                dists[i][j] = helper.get_KL_fast_alt(new_msa, i, j, global_stuff.aa_to_num, weights)
                dists[j][i] = dists[i][j]
        end = datetime.datetime.now()
        print 'calc: ', end - start
        return dists

# right now, only know how to read/write float mats, so make this a regular wrapper
class protein_mutation_list(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['uniprot_id'])

    def whether_to_override(self, object_key):
        return True
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        protein_name = params.get_param('uniprot_id')
        neutral_file = global_stuff.base_folder + protein_name + '/' + 'all_mutations'
        mutations = []
        f = open(neutral_file)

        for line in f:
            if global_stuff.cosmic_or_humvar == 'humvar':
                s = line.strip().split('\t')
            elif global_stuff.cosmic_or_humvar == 'cosmic':
                s = line.strip().split(',')
            try:
                if global_stuff.cosmic_or_humvar == 'humvar':
                    mutations.append([protein_name, int(s[0])-1, s[1], s[2], int(s[3])])
                elif global_stuff.cosmic_or_humvar == 'cosmic':
                    mutations.append([protein_name, int(s[1])-1, s[2], s[3], int(s[4]), int(s[5]), int(s[6])])
            except Exception, err:
                print err


        return mutations

class general_distance(wrapper.mat_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        which_dist = params.get_param('which_dist')
        keys = set(['which_dist'])
        if which_dist == 0:
            return keys | pairwise_dist.get_all_keys(params, self)
        elif which_dist == 1:
            return keys | mf_distance.get_all_keys(params, self)
        elif which_dist == 2:
            return keys | MIP_distance.get_all_keys(params, self)
        elif which_dist == 3:
            return keys | psicov_distance.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):


        which_dist = self.get_param(params,'which_dist')
        if which_dist == 0:
            return self.get_var_or_file(pairwise_dist, params, recalculate, True, False)
        elif which_dist == 1:

            ans = self.get_var_or_file(mf_distance, params, recalculate, True, False)

            return ans
        elif which_dist == 2:
            return self.get_var_or_file(MIP_distance, params)
        elif which_dist == 3:
            return self.get_var_or_file(psicov_distance, params)
        

class mf_distance(wrapper.mat_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        keys = set(['pseudo_c'])
        return keys | general_msa.get_all_keys(params, self) | general_seq_weights.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return False

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        import numpy
        q = global_stuff.q

        msa = self.get_var_or_file(general_msa, params, recalculate, False, False)



        C = numpy.matrix(numpy.zeros(((q-1)*msa.get_alignment_length(),(q-1)*msa.get_alignment_length())),copy=False,dtype=numpy.float32)

        print 'C size: ', C.size


        directed = numpy.zeros((msa.get_alignment_length(), q, msa.get_alignment_length(), q),dtype=numpy.float32)

        print 'directed_size: ', directed.size



        import datetime
        past=datetime.datetime.now()

        print '                   STARTING MF ', params.get_param('uniprot_id')



        


        
        
        # make C and get node counts.  for each edge, need to get joint counts

        weights = self.get_var_or_file(general_seq_weights, params, recalculate, False, False)



        total_weight = numpy.array(weights).sum()

        c = self.get_param(params, 'pseudo_c')

        pseudo_total = c * total_weight
        
        total_weight += pseudo_total

        

        node_counts = numpy.zeros((msa.get_alignment_length(),q)) + (pseudo_total / q)
        node_sizes = numpy.zeros((msa.get_alignment_length())) + (pseudo_total)



        
        mapping = global_stuff.aa_to_num
        for i in range(msa.get_alignment_length()):
            for j in range(len(msa)):
                if msa[j,i] != '-' and msa[j,i] != 'X':
                    node_counts[i,helper.do_map(msa[j,i],mapping)] += weights[j]
                    node_sizes[i] += weights[j]



        for i in range(msa.get_alignment_length()):
            for j in range(q):
                node_counts[i,j] /= node_sizes[i]



        edge_counts = numpy.zeros((msa.get_alignment_length(),q,msa.get_alignment_length(),q),dtype=numpy.uint16) + (pseudo_total / (q*q))
        edge_sizes = numpy.zeros((msa.get_alignment_length(),msa.get_alignment_length()),dtype=numpy.uint16)  + (pseudo_total)
        
        for i in range(msa.get_alignment_length()):
            for j in range(msa.get_alignment_length()):
                for k in range(len(msa)):
                    if msa[k,i] != '-' and msa[k,i] != 'X' and msa[k,j] != '-' and msa[k,j] != 'X':
                        edge_counts[i,helper.do_map(msa[k,i],mapping),j,helper.do_map(msa[k,j],mapping)] += weights[k]
                        edge_sizes[i,j] += weights[k]

        for i in range(msa.get_alignment_length()):
            for j in range(msa.get_alignment_length()):
                for k in range(q):
                    for l in range(q):
                        edge_counts[i,k,j,l] /= edge_sizes[i,j]




        def site_aa_to_index(i, aa):
            return i * (q-1) + aa





        for i in range(msa.get_alignment_length()):
            for j in range(msa.get_alignment_length()):
                for k in range(q-1):
                    for l in range(q-1):
                        try:
                            C[site_aa_to_index(i,k), site_aa_to_index(j,l)] = edge_counts[i,k,j,l] - node_counts[i,k] * node_counts[j,l]
                        except:
                            x=2

                            x=2

        pdb.set_trace()

        C += .0000000005 * numpy.identity(C.shape[0])

        print '                           STARTING INVERSION', datetime.datetime.now() - past
        past = datetime.datetime.now()

        E = C.I
        pdb.set_trace()

        print '                           FINISHED INVERSION', datetime.datetime.now() - past
        past = datetime.datetime.now()




        def get_E(E,i,j,k,l):
            # i and k are positions
            if j == q-1 or l == q-1:
                return 0
            else:
                return -1.0 * E[site_aa_to_index(i,j),site_aa_to_index(k,l)]



        h_node = numpy.zeros((msa.get_alignment_length(), q))

        for i in range(msa.get_alignment_length()):
            mean = 0.0
            for j in range(q):
                temp = math.log(node_counts[i,j])
                for k in range(msa.get_alignment_length()):
                    if i != k:
                        for l in range(q-1):
                            if i == 56 and j == 14 and k == 55:
                                #pdb.set_trace()
                                print '              b ', get_E(E,i,j,k,l), 'a', node_counts[k,l],i,j,k,l
                            temp -= get_E(E,i,j,k,l) * node_counts[k,l]
                h_node[i,j] = temp
                mean += temp




        
        def get_H(H,i,j):
            if j == q-1:
                return 0
            else:
                return H[i,j]



        for i in range(msa.get_alignment_length()):
            for j in range(msa.get_alignment_length()):
                total = 0
                for k in range(q):
                    for l in range(q):
                        try:
                            temp = get_E(E,i,k,j,l) + get_H(h_node,i,k) + get_H(h_node,j,l)
                            if i == 56 and j == 55 and k == 14:
                                print 'c', get_E(E,i,k,j,l) , get_H(h_node,i,k) , get_H(h_node,j,l)
                            directed[i,k,j,l] = temp
                        except:
                            assert False
                            print temp
                        if total == 0:
                            try:
                                total = temp
                            except:
                                pdb.set_trace()
                                x=2
                        else:
                            total = numpy.logaddexp(total, temp)
                asdf = 0


                for k in range(q):
                    for l in range(q):

                        #print directed[i,j,k,l], total
                        if directed[i,k,j,l] > total:
                            #pdb.set_trace()
                            #print directed[i,j,k,l], total, 'bad'
                            pass
                        directed[i,k,j,l] = math.exp(directed[i,k,j,l] - total)


        dist = [ [0.0 for i in range(msa.get_alignment_length())] for j in range(msa.get_alignment_length())]

        # need to get node mar

        for i in range(msa.get_alignment_length()):
            for j in range(msa.get_alignment_length()):
                val = 0
                dir_node_marg_a = numpy.zeros((q,))
                dir_node_marg_b = numpy.zeros((q,))
                for k in range(q):
                    for l in range(q):
                        dir_node_marg_a[k] += directed[i,k,j,l]
                        dir_node_marg_b[l] += directed[i,k,j,l]
                for k in range(q):
                    for l in range(q):
                        if directed[i,k,j,l] > 0:
                            val += directed[i,k,j,l] * (math.log(directed[i,k,j,l]) - math.log(dir_node_marg_a[k]) - math.log(dir_node_marg_b[l]))
                            """
                            try:
                                val += directed[i,j,k,l] * (math.log(directed[i,j,k,l]) - math.log(dir_node_marg[i,k]) - math.log(dir_node_marg[j,l]))
                            except:
                                pdb.set_trace()
                                x=2
                                """

                dist[i][j] = val

        print '                           FINISHED EVERYTHING', datetime.datetime.now() - past
        past = datetime.datetime.now()
        pdb.set_trace()
        return dist
                    
            
                
    

class mutation_list_given_protein_list(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['protein_list_file'])

    def whether_to_override(self, object_key):
        return True
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):



        protein_list_file = self.get_param(params, 'protein_list_file')
        f = open(protein_list_file, 'r')
        mutation_list = []
        i = 0
        for line in f:
            if i % 100 == 0:
                print i
            i += 1

            protein_name = line.strip()
            self.set_param(params, 'uniprot_id', protein_name)

            mutation_list = mutation_list + self.get_var_or_file(protein_mutation_list, params, False, False, False)
        return mutation_list

class general_msa(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        keys = set(['which_msa'])
        which_msa = params.get_param('which_msa')
        if which_msa == 0:
            return keys | agW.get_all_keys(params, self)
        elif which_msa == 1:
            return keys | their_agW.get_all_keys(params, self)
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        which_msa = self.get_param(params, 'which_msa')
        if which_msa == 0:
            return self.get_var_or_file(agW, params, False, False, False)
        elif which_msa == 1:
            return self.get_var_or_file(their_agW, params, False, False, False)

    def whether_to_override(self, object_key):
        return True

class div_weights(wrapper.vect_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return general_msa.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        msa = self.get_var_or_file(general_msa, params, False, False, False)
        return helper.get_weight_of_msa_seqs(msa)

class general_seq_weights(wrapper.vect_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    def whether_to_override(self, object_key):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        keys = set(['which_weight'])
        which_weight = params.get_param('which_weight')
        if which_weight == 1:
            return keys | div_weights.get_all_keys(params, self)
        elif which_weight == 0:
            return keys | general_msa.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        which_weight = self.get_param(params, 'which_weight')
        if which_weight == 1:
            return self.get_var_or_file(div_weights, params, False, False, False)
        elif which_weight == 0:
            msa = self.get_var_or_file(general_msa, params, False, False, False)
            return helper.normalize([1.0 for i in range(len(msa))])
                                        

class edge_to_rank(wrapper.edge_to_int_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return general_distance.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        print "beginning edge_to_rank ", params.get_param('uniprot_id')

        dists = self.get_var_or_file(general_distance, params, recalculate, True, False, False)

        print 'got dists: ', len(dists)

        
        the_dict = {}
        temp = []
        for i in range(len(dists)):
            for j in range(len(dists)):
                if i != j:
                    temp.append([(i,j),dists[i][j]])
                
        temp_sorted = sorted(temp, key = lambda elt: elt[1], reverse = True)
        i = 0
        for elt in temp_sorted:
            the_dict[elt[0]] = i
            i += 1

        print 'finish etr'

        return the_dict


class neighbors_w(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):


        etr = self.get_var_or_file(edge_to_rank, params, recalculate, True, False, False)

        seq = self.get_var_or_file(dW, params, recalculate, True, False, False)

        print 'start to retrieve dists'
        dists = self.get_var_or_file(pairwise_dist, params, recalculate, True, False, False)
        print 'finish to retrieve dists'
        length = len(seq)
        avg_deg = self.get_param(params, 'avg_deg')
        edges = []

        print 'start neighbor'

        rank_cutoff = ((length * 1.0) * avg_deg)
        #pdb.set_trace()
        for i in range(length):
            edges.append([])
            for j in range(length):
                if i != j:
                    if etr[(i,j)] < rank_cutoff:
                        edges[-1].append(j)
                    """
                    if etr[(i,j)] < rank_cutoff:
                        if len(edges[-1]) == 0:
                            edges[-1].append(j)
                        else:
                            if etr[(i,j)] < etr[(i,edges[-1][0])]:
                                edges[-1][0] = j
                    """
        #pdb.set_trace()

        print 'end neighbor'

        
        return edges

class neighbors_w_weight_w(wrapper.int_float_tuple_mat_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        keys = set(['avg_deg'])
        return keys | edge_to_rank.get_all_keys(params, self) | dW.get_all_keys(params, self) | general_msa.get_all_keys(params, self) | dW.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return False

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):


        etr = self.get_var_or_file(edge_to_rank, params, recalculate, True, False, False)

        seq = self.get_var_or_file(dW, params, recalculate, True, False, False)

        print 'start to retrieve dists'
        dists = self.get_var_or_file(general_distance, params, recalculate, True, False, False)
        print 'finish to retrieve dists'
        length = len(seq)
        avg_deg = self.get_param(params, 'avg_deg')
        edges = []

        print 'start neighbor'

        rank_cutoff = ((length * 1.0) * avg_deg)

        def frac_non_skip(msa):
            ans = []
            for i in range(msa.get_alignment_length()):
                ans.append((len(msa) - msa.get_column(i).count('-'))/float(len(msa)))
            return ans

        def ok(pct, i, j):
            if (pct[i] + pct[j])/2.0 < 0:
                return False
            else:
                return True

        # figure out how many edges to throw out

        msa = self.get_var_or_file(general_msa, params, False, False, False)
        pct = frac_non_skip(msa)
        num_bad = 0
        for i in range(length):
            for j in range(length):
                try:
                    if not ok(pct,i,j):
                        num_bad += 1
                except:
                    x=2
                    pdb.set_trace()

        rank_cutoff += num_bad

        for i in range(length):
            edges.append([])
            for j in range(length):
                if i != j:
                    try:
                        if etr[(i,j)] < rank_cutoff and ok(pct,i,j):
                            edges[-1].append((j,dists[i][j]))
                    except:
                        x=2
                        pdb.set_trace()
                    """
                    if etr[(i,j)] < rank_cutoff:
                        if len(edges[-1]) == 0:
                            edges[-1].append(j)
                        else:
                            if etr[(i,j)] < etr[(i,edges[-1][0])]:
                                edges[-1][0] = j
                    """
        #pdb.set_trace()

        print 'end neighbor'

        
        return edges






class filtered_mutation_list_given_protein_list(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return mutation_list_given_protein_list.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return False

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        all_mutations = self.get_var_or_file(mutation_list_given_protein_list, params, recalculate, False, False)

        num = 0
        num_add = 0

        n_dict = {}
        q_dict = {}
        m_dict = {}

        neighbor_cutoff = self.get_param(params, 'n_cutoff')
        filter_cutoff = self.get_param(params, 'f_cutoff')

        filtered_mutations = []

        for mutation in all_mutations:

            protein_name = mutation[0]
            #print protein_name
            self.set_param(params, 'uniprot_id', protein_name)

            seq = self.get_var_or_file(dW, params, recalculate, False, False)

            """
                try:
                    neighbors = n_dict[protein_name]
                except:
                    neighbors = self.get_var_or_file(neighbors_w, params, recalculate, True, False)
                    n_dict[protein_name] = neighbors

                try:
                    query = q_dict[protein_name]
                except:
                    query = self.get_var_or_file(dW, params, recalculate, True, False)
                    q_dict[protein_name ] = query
                    """
            """
                try:
                    neighbors = n_dict[protein_name]
                except:
                    neighbors = self.get_var_or_file(neighbors_w_weight_w, params, recalculate, True, False)
                    n_dict[protein_name] = neighbors
                
                try:
                    msa = m_dict[protein_name]
                except:
                    msa = self.get_var_or_file(general_msa, params, recalculate, True, False)
                    m_dict[protein_name ] = msa

                    """
            if len(seq) < 500000:
                neighbors = self.get_var_or_file(neighbors_w_weight_w, params, recalculate, True, False)
                msa = self.get_var_or_file(general_msa, params, recalculate, True, False)

                try:    
                    pos = mutation[1]


                    col = msa.get_column(pos)
                
                    if col.count(mutation[3]) > 5 and len(msa) > 200 and len(neighbors[pos])>0:
                    #if True:
                        #print 'added: ', protein_name
                        filtered_mutations.append(mutation)
                        num_add += 1
                except Exception, err:
                    print err
                num += 1
                if num % 50 == 0:
                    print '                               ', num_add, num

                """
                filtered_column = helper.filter_msa_based_on_pos_neighbors_and_query(query, pos, msa, neighbors[pos])

                print '       ', len(neighbors[pos]), len(filtered_column), neighbor_cutoff, filter_cutoff

                if len(neighbors[pos]) > neighbor_cutoff and len(filtered_column) > filter_cutoff:
                    filtered_mutations.append(mutation)
                    print '                                   YES'
                """


                
        return filtered_mutations

def get_protein_info(protein_list, info_file, params):
    import wc
    f = open(protein_list, 'r')
    g = open(info_file, 'w')
    i = 0
    for line in f:
        name = line.strip()
        print name,i
        i += 1
        try:
            params.set_param('uniprot_id', name)

            def get_name():
                return [name]

            def get_overlap():
                params.set_param('which_dist', 0)
                n1 = wc.get_stuff(general_distance, params, False, False, False)
                params.set_param('which_dist', 1)
                n2 = wc.get_stuff(general_distance, params, False, False, False)
                n1set = set()
                n2set = set()
                for i in len(seq):
                    for x in n1[i]:
                        n1set.add((i,x[0]))
                    for x in n2[i]:
                        n2set.add((i,x[1]))
                both = n1set.intersection(n2set)
                frac = len(both) / float(len(n1set))
                return [str(both), str(frac)]

            def get_seq_length():
                seq = wc.get_stuff(dW, params, False, False, False)
                return [str(len(seq))]
            
            def get_msa_length():
                msa = wc.get_stuff(general_msa, params, False, False, False)
                return [str(len(msa))]

            def get_num_mutations():
                mutations = wc.get_stuff(protein_mutation_list, params, False, False, False)
                return [str(len(mutations))]

            info = []
            #which_info = [get_name, get_seq_length, get_msa_length, get_num_mutations]
            which_info = [get_name, get_seq_length, get_msa_length]
            for which in which_info:
                info = info + which()
            
            g.write(string.join(info,sep=',') + '\n')
            
        except:

            print 'fail'
    f.close()
    g.close()


def get_every_site_info(params, protein_list, dist_file):
    import wc
    f = open(dist_file,'w')

    protein_list = protein_list
    h = open(protein_list)

    k=0

    for line in h:
        protein_name = line.strip()
        print protein_name, k
        k += 1
        params.set_param('uniprot_id', protein_name)
        #params.set_param('which_dist', 0)
        #pdb.set_trace()
        params.set_param('which_dist',1)
        all_neighbors = {}
        avg_degs = range(1,6)
        for i in avg_degs:
            params.set_param('avg_deg',i)
            all_neighbors[i] = wc.get_stuff(neighbors_w_weight_w, params, False, False, False)
        

        msa = wc.get_stuff(their_agW, params, False, False, False)
        for i in range(len(all_neighbors[avg_degs[0]])):
            not_skip = len(msa) - msa.get_column(i).count('-')
            f.write(str(i) + ',' + protein_name + ',' + str(len(msa)) + ',' + str(not_skip) + ',' + str(float(not_skip)/len(msa)))
            for j in avg_degs:

                f.write(','+str(len(all_neighbors[j][i])))
            f.write('\n')


    """

        msa = wc.get_stuff(their_agW, params, False, False, False)
        dists = wc.get_stuff(general_distance, params, False, False, False)
        no_skips = []
        for i in range(msa.get_alignment_length()):
            col = msa.get_column(i)
            no_skips.append([x for x in col if x != '-'])

        


    temp = []
    for i in range(len(neighbors)):
        for it in neighbors[i]:
            pdb.set_trace()
            j = it[0]
            temp.append([i,j,it[1],len(no_skips[i]),len(no_skips[j])])

    sorted_temp = sorted(temp, key = lambda x:x[2])
    g = open(len_file)
    import string
    for x in sorted_temp:
        g.write(string.join([str(it) for it in x],',') + '\n')
    g.close()

    """

# for computations involving mutations, specify protein_list to draw from.

def get_overlap(params, protein_list, out_file):


    f = open(protein_list, 'r')
    import wc
    for line in f:
        name = line.strip()
        params.set_param('uniprot_id', name)

        params.set_param('which_dist', 0)
        n1 = wc.get_stuff(neighbors_w_weight_w, params, False, False, False)
        params.set_param('which_dist', 1)
        n2 = wc.get_stuff(neighbors_w_weight_w, params, False, False, False)

        print helper.get_overlap(n1,n2)

    


def get_mutation_info(protein_list_file, out_file, params):
    import wc

    params.set_param('protein_list_file', protein_list_file)
    pdb.set_trace()
    l = wc.get_stuff(mutation_list_given_protein_list, params, False, False, False, False)

    f = open(out_file, 'w')

    avg_degs = [2]
    
    i = 0
    for mutation in l:

        name = mutation[0]
        params.set_param('uniprot_id', name)

        def get_name():
            return [mutation[0]]

        def get_pos():
            return [str(mutation[1])]

        def get_deg():
            ans = []
            for deg in avg_degs:
                params.set_param('avg_deg', deg)
                all_neighbors = wc.get_stuff(neighbors_w_weight_w, params, False, False, False)
                ans.append(str(len(all_neighbors[mutation[1]])))
            return ans

        def get_wild_num():
            msa = wc.get_stuff(general_msa, params, False, False, False)
            col = msa.get_column(mutation[1])
            return [str(col.count(mutation[2]))]

        def get_wild_category_num():
            msa = wc.get_stuff(general_msa, params, False, False, False)
            col = msa.get_column(mutation[1])
            wild_res_cat = global_stuff.aa_to_class[mutation[2]]
            count = 0
            for aa in col:
                if global_stuff.aa_to_class[aa] == wild_res_cat:
                    count += 1
            return [str(count)]

        def get_mut_category_num():
            msa = wc.get_stuff(general_msa, params, False, False, False)
            col = msa.get_column(mutation[1])
            mut_res_cat = global_stuff.aa_to_class[mutation[3]]
            count = 0
            for aa in col:
                if global_stuff.aa_to_class[aa] == mut_res_cat:
                    count += 1
            return [str(count)]

        def get_is_same_cat():
            if global_stuff.aa_to_class[mutation[3]] == global_stuff.aa_to_class[mutation[2]]:
                return ['1']
            else:
                return ['0']
        
        def get_mut_num():
            msa = wc.get_stuff(general_msa, params, False, False, False)
            col = msa.get_column(mutation[1])
            return [str(col.count(mutation[3]))]

        def get_msa_length():
            msa = wc.get_stuff(general_msa, params, False, False, False)
            return [str(len(msa))]

        def get_cosmic_info():
            return [str(mutation[4]), str(mutation[5]), str(mutation[6])]

        def get_whether_bad():
            return [str(mutation[4])]

        info = []
        if i % 50 == 0:
            print get_name(), i
        i += 1
        if global_stuff.cosmic_or_humvar == 'cosmic':
            which_info = [get_name, get_pos, get_wild_num, get_mut_num, get_deg, get_cosmic_info, get_msa_length]
        elif global_stuff.cosmic_or_humvar == 'humvar':
            which_info = [get_name, get_pos, get_wild_num, get_mut_num, get_deg, get_whether_bad, get_msa_length]
        try:
            for which in which_info:
                info = info + which()
            
            f.write(string.join(info,sep=',') + '\n')
        except Exception, err:
            print err
            pass
        

    f.close()

# outputs file.  could be roc file input, or other input.  one argument is function that assigns a number to each mutation.  the function determines the output
def get_output_file(params, protein_list_file, out_file, use_neighbor, ignore_pos, max_neighbors, num_trials, pseudo_total, sim_f, norm_f, mut_to_num_f):
    import wc
    params.set_param('protein_list_file', protein_list_file)
    l = wc.get_stuff(filtered_mutation_list_given_protein_list, params)
    i = 0
    scores = []
    labels = []
    bad_count = 0

    for mutation in l:
        if i%1 == 0:
            print i, 'calculating'

        i += 1
        try:

            score = helper.predict_position_energy_weighted(params, mutation, use_neighbor, ignore_pos, max_neighbors, num_trials, pseudo_total, sim_f)
            scores.append(score)
            labels.append(mut_to_num_f(mutation))
        except Exception, err:
            bad_count += 1
            pdb.set_trace()
            print err, mutation, bad_count



    assert len(scores) == len(labels)

    normed_scores = norm_f(scores)
    ans = [ [0, 0] for i in range(len(scores))]
    for i in range(len(scores)):
        ans[i][0] = normed_scores[i]
        ans[i][1] = labels[i]

    helper.write_mat(ans, out_file)


    
    

def get_roc_file(params, in_file, out_file, use_neighbor, ignore_pos, max_neighbors, weighted, num_trials, pseudo_total, sim_f):
    import wc
    recalculate = False


    params.set_param('protein_list_file', in_file)
    
    l = wc.get_stuff(filtered_mutation_list_given_protein_list, params, recalculate, False, False, False)

    ans = []
    i = 0
    bad_count = 0
    for mutation in l:
        if i % 100 == 0:
            print '      ', i


        i += 1
        try:
            score = helper.predict_position_energy_weighted(params, recalculate, mutation, use_neighbor, ignore_pos, max_neighbors, weighted, num_trials, pseudo_total, sim_f)
            ans.append([score, mutation[4]])
        except Exception, err:
            bad_count += 1
            print 'failed to get score', err, bad_count

    helper.write_mat(ans, out_file)



