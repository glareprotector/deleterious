
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
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

import math
import subprocess
import string
import os
import random
import pdb
import re


import global_stuff

import sys

import wrapper

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
        return True

# blast results file wrapper(xml format)
class adW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        if params.get_param('which_blast') == 0:
            return set(['ev', 'which_blast', 'psiblast_iter']) | bW.get_all_keys(params, self)
        else:
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

        print >> sys.stderr, 'RUNNING BLAST!!!!!!!'
        if self.get_param(params, 'which_blast') == 0:

            psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLAST_PATH, outfmt = 5, query = '\''+f.name+'\'', db = global_stuff.BLASTDB_PATH, out = self.get_holding_location(), evalue = self.get_param(params, 'ev'), num_iterations = self.get_param(params, 'psiblast_iter'))

        elif self.get_param(params, 'which_blast') == 1:
            #psi_blast_cline = global_stuff.BLASTP_PATH + 
            psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLASTP_PATH, outfmt = 5, query = '\''+f.name+'\'', db = global_stuff.BLASTDB_PATH, out = self.get_holding_location(), evalue = self.get_param(params, 'ev'))
        elif self.get_param(params, 'which_blast') == 2:
            psi_blast_cline = global_stuff.DELTABLAST_PATH + ' -outfmt ' + '5' + ' -evalue ' + str(self.get_param(params, 'ev')) + ' -query ' + '\''+f.name+'\'' + ' -db ' + global_stuff.BLASTDB_PATH + ' -out ' + self.get_holding_location() + ' -rpsdb ' + global_stuff.CDD_PATH
        #psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLAST_PATH, outfmt = 5, query = '\''+f.name+'\'', out = self.get_holding_location())
        #pdb.set_trace()
        print >> sys.stderr, psi_blast_cline
        pdb.set_trace()
        subprocess.call(str(psi_blast_cline), shell=True, stdout = sys.stderr, stderr = sys.stderr, executable='/bin/bash')

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
        try:
            record = NCBIXML.read(open(blast_xml_handle.name, 'r'))
        except Exception, err:
            records = NCBIXML.parse(open(blast_xml_handle.name, 'r'))
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
        print >> sys.stderr, 'WROTE ', self.get_holding_location()

        return output_handle


# gets the result of msa in fasta file.  gaps not removed yet.  query should have name of 'QUERY'
class afW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    def whether_to_override(self, object_key):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        if params.get_param('which_msa') == 0:
            return set(['which_msa']) | aeW.get_all_keys(params, self)
        elif params.get_param('which_msa') == 2:
            return set(['which_msa']) | hhblits_msa_file.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        if self.get_param(params, 'which_msa') == 0:
            msa_input_handle = self.get_var_or_file(aeW, params, recalculate, to_pickle, to_filelize, always_recalculate)
            cline = MuscleCommandline(cmd = global_stuff.MUSCLE_PATH, input = '\''+msa_input_handle.name+'\'', out = self.get_holding_location(), clw = False, maxiters = 2)
            subprocess.call(str(cline), shell=True, stdout = sys.stderr, stderr = sys.stderr, executable='/bin/bash')
            return open(self.get_holding_location())
        elif self.get_param(params, 'which_msa') == 2:

            #get the a3m file, convert to fasta, read in, find query, write out msa but with query renamed
            hhblits_msa = self.get_var_or_file(hhblits_msa_file, params)
            temp_fasta_f = self.get_holding_location() + '.temp_fasta'
            convert_cmd = global_stuff.HHBLITS_CONVERT_A3M_TO_FASTA + ' a3m ' + ' fas ' + '\''+hhblits_msa.name+'\'' + ' ' + temp_fasta_f
         
            try:
                code = subprocess.call(convert_cmd, shell=True, stdout = sys.stderr, stderr = sys.stderr, executable='/bin/bash')
            except Exception, err:
                print err


            print code
            
            temp_msa = AlignIO.read(temp_fasta_f, 'fasta')

            try:
                import os
                os.remove(temp_fasta_f)
            except Exception, err:
                print err

            query_seq = self.get_var_or_file(dW, params)
            query_name = query_seq.name
            

            sequences_to_write = []
            
            for seq in temp_msa:
                if seq.id != query_name:
                    sequences_to_write.append(seq)
                else:
                    sequences_to_write.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq.seq.tostring()), id = 'QUERY'))

            output_handle = open(self.get_holding_location(), 'w')
            SeqIO.write(sequences_to_write, output_handle, 'fasta')
            return output_handle


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
        print cmd

        subprocess.call(cmd, stdout = sys.stderr, stderr = sys.stderr, shell=True, executable='/bin/bash')
        return open(self.get_holding_location(), 'r')

    def whether_to_override(self, object_key):
        return False

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
        print >> sys.stderr, cmd
        subprocess.call(cmd, stdout = sys.stderr, stderr = sys.stderr, shell=True, executable='/bin/bash')
        subprocess.call('mv ' + self.get_holding_location() + '_MIp.txt ' + self.get_holding_location(), stdout = sys.stderr, stderr = sys.stderr, shell=True, executable='/bin/bash')
        subprocess.call('rm ' + self.get_holding_location() + '_count.txt', shell=True, stdout = sys.stderr, stderr = sys.stderr, executable='/bin/bash')
        subprocess.call('rm ' + self.get_holding_location() + '_MIp.dot', shell=True, stdout = sys.stderr, stderr = sys.stderr, executable='/bin/bash')
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
        return general_afW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        print 'AGW!!!!!'

        msa = self.get_var_or_file(general_afW, params, recalculate, False, False, False)
        # msa = AlignIO.read(f.name, 'fasta')
        # search for the query sequence
        idx = -1
        for i in range(len(msa)):
            if msa[i].id == 'QUERY':
                idx = i
                break


        seq = self.get_var_or_file(dW, params)

        seqs = [ [None for i in range(len(seq))] for j in range(len(msa)) ]


        pos = 0
        for i in range(msa.get_alignment_length()):
            if msa[idx,i] != '-':
               for j in range(len(msa)):

                   seqs[j][pos] = msa[j,i]
               pos += 1


        to_add = [ SeqRecord(Seq(string.join(seqs[i],sep='')), id = msa[i].id) for i in range(len(msa))]
        return MultipleSeqAlignment(to_add)

            
        # find the first non-insertion column
        i = 0

        while msa[idx,i] == '-':
            #print >> sys.stderr, msa[idx,i]
            i = i + 1
            #pdb.set_trace()

        to_return = msa[:,i:(i+1)]
        # add in all the other columns

        for k in range(i+1, msa.get_alignment_length()):
            print k
            if msa[idx,k] != '-':
                #print >> sys.stderr, k
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
            #print >> sys.stderr, len(to_write), to_write
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
        print >> sys.stderr, 'length: ', msa.get_alignment_length(), len(msa)
        import pdb
        import datetime
        start = datetime.datetime.now()
        dists = [ [0 for j in range(msa.get_alignment_length())] for i in range(msa.get_alignment_length())]

        new_msa = [ [msa[i][j] for j in range(msa.get_alignment_length())] for i in range(len(msa))]

        weights = self.get_var_or_file(general_seq_weights, params, False, False, False)



        print >> sys.stderr, datetime.datetime.now() - start

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
        print >> sys.stderr, 'calc: ', end - start
        return dists

# right now, only know how to read/write float mats, so make this a regular wrapper
class humvar_protein_mutation_list(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['uniprot_id'])

    def whether_to_override(self, object_key):
        return False
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        protein_name = params.get_param('uniprot_id')
        neutral_file = global_stuff.base_folder + protein_name + '/' + 'all_mutations'
        mutations = []
        f = open(neutral_file)

        for line in f:
            s = line.strip().split('\t')
            try:
                mutations.append([protein_name, int(s[0])-1, s[1], s[2], int(s[3])])
            except Exception, err:
                print >> sys.stderr, err


        return mutations


class saapdb_protein_mutation_list(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

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
            s = line.strip().split('\t')
            try:
                mutations.append([protein_name, int(s[1])-1, s[2], s[3], int(s[4])])
            except Exception, err:
                print >> sys.stderr, err


        return mutations


class cosmic_protein_mutation_list(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['uniprot_id'])

    def whether_to_override(self, object_key):
        return False
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        protein_name = params.get_param('uniprot_id')
        neutral_file = global_stuff.base_folder + protein_name + '/' + 'all_mutations'
        mutations = []
        f = open(neutral_file)

        for line in f:
            s = line.strip().split(',')
            try:
                mutations.append([protein_name, int(s[1])-1, s[2], s[3], int(s[4]), int(s[5]), int(s[6])])
            except Exception, err:
                print >> sys.stderr, err


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
        return False

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

        print >> sys.stderr, 'C size: ', C.size


        directed = numpy.zeros((msa.get_alignment_length(), q, msa.get_alignment_length(), q),dtype=numpy.float32)

        print >> sys.stderr, 'directed_size: ', directed.size



        import datetime
        past=datetime.datetime.now()

        print >> sys.stderr, '                   STARTING MF ', params.get_param('uniprot_id')



        


        
        
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

        print >> sys.stderr, '                           STARTING INVERSION', datetime.datetime.now() - past
        past = datetime.datetime.now()

        E = C.I
        pdb.set_trace()

        print >> sys.stderr, '                           FINISHED INVERSION', datetime.datetime.now() - past
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
                                print >> sys.stderr, '              b ', get_E(E,i,j,k,l), 'a', node_counts[k,l],i,j,k,l
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
                                print >> sys.stderr, 'c', get_E(E,i,k,j,l) , get_H(h_node,i,k) , get_H(h_node,j,l)
                            directed[i,k,j,l] = temp
                        except:
                            assert False
                            print >> sys.stderr, temp
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

                        #print >> sys.stderr, directed[i,j,k,l], total
                        if directed[i,k,j,l] > total:
                            #pdb.set_trace()
                            #print >> sys.stderr, directed[i,j,k,l], total, 'bad'
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

        print >> sys.stderr, '                           FINISHED EVERYTHING', datetime.datetime.now() - past
        past = datetime.datetime.now()
        pdb.set_trace()
        return dist
                    
            
class humvar_mutation_list_all(wrapper.obj_wrapper):

    def whether_to_override(self, object_key):
        return True
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(global_stuff.humvar_genes_file, 'r')
        mutation_list = []
        for line in f:
            protein_name = line.strip()
            self.set_param(params, 'uniprot_id', protein_name)
            mutation_list = mutation_list + self.get_var_or_file(humvar_protein_mutation_list, params)
        return mutation_list

class current_cosmic_mutation_list_all(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()

    def whether_to_override(self, object_key):
        return True
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        # read in a list of cosmic gene names 
        f = open(global_stuff.cosmic_genes_file, 'r')
        mutation_list = []
        i = 0
        for line in f:
            if i % 100 == 0:
                print >> sys.stderr, i
            i += 1

            protein_name = line.strip()
            self.set_param(params, 'uniprot_id', protein_name)

            mutation_list = mutation_list + self.get_var_or_file(cosmic_protein_mutation_list, params, False, False, False)
        return mutation_list




class saapdb_mutation_list_all(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()

    def whether_to_override(self, object_key):
        return True
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        # read in a list of cosmic gene names 
        f = open(global_stuff.saapdb_genes_file, 'r')
        mutation_list = []
        i = 0
        for line in f:
            if i % 100 == 0:
                print >> sys.stderr, i
            i += 1

            protein_name = line.strip()
            self.set_param(params, 'uniprot_id', protein_name)

            mutation_list = mutation_list + self.get_var_or_file(saapdb_protein_mutation_list, params, False, False, False)
        return mutation_list


class p53_mutation_list_all(wrapper.obj_wrapper):

    def whether_to_override(self, object_key):
        return True

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(global_stuff.P53_MUTATIONS, 'r')
        f.next()
        mutation_list = []
        pdb.set_trace()
        for line in f:
            s = line.strip().split('\t')
            raw = s[1]
            wild = raw[0]
            mut = raw[-1]
            pos = int(raw[1:-1])
            mutation = ['TP53', pos-1, wild, mut, float(s[3]), float(s[4]), float(s[5]), float(s[6]), float(s[7]), float(s[8]), float(s[9]), float(s[10])]
            mutation_list.append(mutation)
        return mutation_list


class reverse_stone_mutation_list_all(wrapper.obj_wrapper):

    def whether_to_override(self, object_key):
        return True

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()


    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):



        params.set_param('uniprot_id', global_stuff.REVERSE_GENE_NAME_TO_USE)


        
        # to go from their file to position in their sequence, subtract 43
        # 2 files: their file, their seq, and the seq i will use to get neighbors
        short_seq = SeqIO.read(global_stuff.REVERSE_SHORT_FILE,'fasta')
        long_seq = self.get_var_or_file(dW, params)
        f = open(global_stuff.REVERSE_MUT_FILE, 'r')
        f.next()

        #align short and long sequence
        short_to_long, long_to_short = helper.get_alignment_mapping(short_seq.seq.tostring(), long_seq.seq.tostring())
        mutation_list = []
        for line in f:
            s = line.strip().split('\t')
            wild = s[0]
            short_pos = int(s[1]) - 44
            mut = s[2]
            cat = line.count('+')
            try:
                long_pos = short_to_long[short_pos]
                assert long_seq[long_pos] == wild
                mutation = [global_stuff.REVERSE_GENE_NAME_TO_USE, long_pos, wild, mut, cat]
                mutation_list.append(mutation)
            except:
                pass
        return mutation_list

class lys_stone_mutation_list_all(wrapper.obj_wrapper):

    def whether_to_override(self, object_key):
        return True

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(global_stuff.LYSOZYME_FILE, 'r')
        f.next()
        mutation_list = []
        for line in f:
            s = line.strip().split()
            pos = int(s[1]) - 1
            wild = s[0]
            mut = s[2]
            cls = line.count('+')
            mutation = ['LYSOZYME', pos, wild, mut, cls]
            mutation_list.append(mutation)
        return mutation_list

class hhblits_msa_file(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    def whether_to_override(self, object_key):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['hhblits_iter', 'ev']) | bW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        f = self.get_var_or_file(bW, params, False, False, False)

        cmd = global_stuff.HHBLITS_PATH + ' -i ' + '\''+f.name+'\'' + ' -d ' + global_stuff.HHBLITS_DB_PATH + ' -oa3m ' + self.get_holding_location() + ' -cpu ' + str(1) + ' -n ' + str(self.get_param(params, 'hhblits_iter')) + ' -e ' + str(self.get_param(params, 'ev'))
        subprocess.call(cmd, shell=True, stderr = sys.stderr, stdout = sys.stderr, executable='/bin/bash')

        return open(self.get_holding_location(), 'r')


        
class general_msa(wrapper.msa_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        keys = set(['which_msa'])
        which_msa = params.get_param('which_msa')
        if which_msa == 1:
            return keys | their_agW.get_all_keys(params, self)
        else:
            return keys | agW.get_all_keys(params, self)
        
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        which_msa = self.get_param(params, 'which_msa')
        if which_msa == 1:
            return self.get_var_or_file(their_agW, params, False, False, False)
        else:
            return self.get_var_or_file(agW, params, False, False, False)
        


    def whether_to_override(self, object_key):
        return True

class fake_cluster_file(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return renamed_afW.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(self.get_holding_location(), 'w')
        msa = self.get_var_or_file(renamed_general_msa, params)
        num_seq = len(msa)
        f.write('Number of clusters : ' + str(0) + '\n')
        f.write('\n')
        f.write('Cluster 0 ; size=' + str(num_seq) + '\n')
        for i in range(num_seq):
            f.write(msa[i].id[0:30] + '\n')
        f.close()
        return open(self.get_holding_location())

class renamed_afW(wrapper.msa_obj_wrapper, wrapper.by_uniprot_id_wrapper):
    """
    renames the sequence id's with numbers
    """

    @classmethod
    def get_all_keys(cls, params, self=None):
        return afW.get_all_keys(params, self)
    
    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        f = self.get_var_or_file(afW, params)
        msa = AlignIO.read(f.name, 'fasta')
        count = 0
        for seq in msa:
            if seq.id != 'QUERY':
                seq.id = (str(count) + '_' + seq.id)[0:30]
                seq.name = seq.id
                count += 1

        return msa

        l = len(msa)

        query = msa[0]
        ans = msa[500:]
        ans.append(query)
        return ans

    def whether_to_override(self, params):
        return True


class leoned_general_msa(wrapper.msa_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    def whether_to_override(self, object_key):
        return True

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['to_cluster']) | renamed_general_msa.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        msa = self.get_var_or_file(renamed_general_msa, params)
        import wc
        msa_file_name = wc.get_wrapper_instance(renamed_general_msa).get_file_location(params)

        working_file = self.get_holding_location() + '_msa'
        query_name = 'QUERY'
        log_file = '/dev/null'
        to_cluster = self.get_param(params, 'to_cluster')

        subprocess.call('cp' +  ' ' + '\''+msa_file_name+'\'' + ' ' + working_file, stdout = sys.stderr, stderr = sys.stderr, shell=True, executable='/bin/bash')
        if to_cluster == 1:
            subprocess.call(string.join([global_stuff.LEON_PATH, working_file, query_name, self.get_holding_location(), working_file, log_file, str(1)],sep=' '), stdout = sys.stderr, stderr = sys.stderr, shell=True, executable = '/bin/bash')
        elif to_cluster == 0:

            cluster_file_name = wc.get_wrapper_instance(fake_cluster_file).get_file_location(params)
            temp_cluster_location = working_file + '.clu'
            subprocess.call('cp' + ' ' + '\''+cluster_file_name+'\'' + ' ' + temp_cluster_location, shell=True, executable='/bin/bash')
            subprocess.call(string.join([global_stuff.LEON_PATH, working_file, query_name, self.get_holding_location(), working_file, log_file, str(0)],sep=' '), stdout = sys.stderr, stderr = sys.stderr, shell=True, executable = '/bin/bash')
        return AlignIO.read(self.get_holding_location(), 'fasta')


class rascalled_afW(wrapper.msa_obj_wrapper, wrapper.by_uniprot_id_wrapper):


    def whether_to_override(self, object_key):
        return True

    @classmethod
    def get_all_keys(cls, params, self=None):
        return renamed_afW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        msa = self.get_var_or_file(renamed_afW, params)
        import wc
        print 'RASCAL!!!'
        msa_file_name = wc.get_wrapper_instance(renamed_afW).get_file_location(params)
        working_file = self.get_holding_location() + '_msa'
        subprocess.call('cp' +  ' ' + '\''+msa_file_name+'\'' + ' ' + working_file, stdout = sys.stderr, stderr = sys.stderr, shell=True, executable='/bin/bash')
        subprocess.call(global_stuff.RASCAL_PATH + ' ' + working_file + ' ' + self.get_holding_location(), stdout = sys.stderr, stderr = sys.stderr, shell=True, executable='/bin/bash')
        return AlignIO.read(self.get_holding_location(), 'fasta')


class norMD_afW(wrapper.msa_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    def whether_to_override(self, object_key):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        to_rascal = params.get_param('to_rascal')

        if to_rascal == 0:
            return set(['to_rascal', 'norm_co']) | rascalled_afW.get_all_keys(params, self)
        elif to_rascal == 1:
            return set(['to_rascal', 'norm_co']) | renamed_afW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        to_rascal = self.get_param(params, 'to_rascal')
        print 'NORMD!!!!!'
        import wc
        if to_rascal == 1:
            msa = self.get_var_or_file(rascalled_afW, params)
            msa_file_name = wc.get_wrapper_instance(rascalled_afW).get_file_location(params)
        elif to_rascal == 0:
            msa = self.get_var_or_file(renamed_afW, params)
            msa_file_name = wc.get_wrapper_instance(renamed_afW).get_file_location(params)
        in_msa_file = self.get_holding_location() + '_msa'
        helper.cp_file(msa_file_name, in_msa_file)
        in_msa_gcg_file =  in_msa_file + '_gcg'
        helper.conv_seq(in_msa_file, in_msa_gcg_file, 'gcg')
        out_msa_gcg_file = self.get_holding_location() + '_out_gcg'
        norm_co = self.get_param(params, 'norm_co')
        subprocess.call(global_stuff.NORMD_PATH + ' ' + in_msa_gcg_file + ' ' + global_stuff.NORMD_SIMMAT + ' ' + str(12) + ' ' + str(1) + ' ' + str(norm_co) + ' ' + out_msa_gcg_file + ' QUERY', stdout = sys.stderr, stderr = sys.stderr, shell=True, executable='/bin/bash')
        #subprocess.call(global_stuff.NORMD_PATH + ' ' + in_msa_gcg_file + ' ' + global_stuff.NORMD_SIMMAT + ' ' + str(12) + ' ' + str(1) + ' ' + str(9.0) + ' ' + 'QUERY' + ' ' + out_msa_gcg_file, shell=True, stdout = sys.stderr, stderr = sys.stderr, executable='/bin/bash')

        helper.conv_seq(out_msa_gcg_file, self.get_holding_location(), 'fasta')
        ans = AlignIO.read(self.get_holding_location(), 'fasta')
        subprocess.call('rm ' + self.get_holding_location() + '*', stdout = sys.stderr, stderr = sys.stderr, shell=True, executable='/bin/bash')
        #helper.rm_files([in_msa_file, in_msa_gcg_file, out_msa_gcg_file, self.get_holding_location()])
        return ans
        
class general_afW(wrapper.msa_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):


        to_rascal = params.get_param('to_rascal')
        to_normd = params.get_param('to_normd')
        if to_rascal == 0 and to_normd == 0:
            return set(['to_rascal','to_normd']) | renamed_afW.get_all_keys(params, self)
        elif to_normd == 0 and to_rascal == 1:
            return set(['to_rascal','to_normd']) | rascalled_afW.get_all_keys(params, self)
        elif to_normd == 1:
            return set(['to_rascal','to_normd']) | norMD_afW.get_all_keys(params, self)
        else:
            raise Exception

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        to_rascal = self.get_param(params, 'to_rascal')
        to_normd = self.get_param(params, 'to_normd')
        if to_rascal == 0 and to_normd == 0:
            return self.get_var_or_file(renamed_afW, params)
        elif to_normd == 0 and to_rascal == 1:
            return self.get_var_or_file(rascalled_afW, params)
        elif to_normd == 1:
            return self.get_var_or_file(norMD_afW, params)
        else:
            raise Exception
    

class div_weights(wrapper.vect_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return general_msa.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        print 'WEIGHTS!!!'
        import wrapper
        msa = self.get_var_or_file(wrapper.my_msa_obj_wrapper, params, False, False, False)
        return helper.get_weight_of_msa_seqs(msa)

class general_seq_weights(wrapper.vect_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    def whether_to_override(self, object_key):
        return True

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

        print >> sys.stderr, "beginning edge_to_rank ", params.get_param('uniprot_id')

        dists = self.get_var_or_file(general_distance, params, recalculate, True, False, False)

        print >> sys.stderr, 'got dists: ', len(dists)

        
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

        print >> sys.stderr, 'finish etr'

        return the_dict


class neighbors_w(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):


        etr = self.get_var_or_file(edge_to_rank, params, recalculate, True, False, False)

        seq = self.get_var_or_file(dW, params, recalculate, True, False, False)

        print >> sys.stderr, 'start to retrieve dists'
        dists = self.get_var_or_file(pairwise_dist, params, recalculate, True, False, False)
        print >> sys.stderr, 'finish to retrieve dists'
        length = len(seq)
        avg_deg = self.get_param(params, 'avg_deg')
        edges = []

        print >> sys.stderr, 'start neighbor'

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

        print >> sys.stderr, 'end neighbor'

        
        return edges

class general_neighbors_w_weight_w(wrapper.int_float_tuple_mat_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    def whether_to_override(self, object_key):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        which_neighbor = params.get_param('which_neighbors')
        if which_neighbor == 0:
            return (set(['which_neighbors']) | neighbors_w_weight_w.get_all_keys(params, self)) 
        elif which_neighbor == 1:
            return (set(['which_neighbors']) | start_neighbors_from_pdb.get_all_keys(params, self)) 
                     


    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        which_neighbors = self.get_param(params, 'which_neighbors')

        if which_neighbors == 0:
            return self.get_var_or_file(neighbors_w_weight_w, params)
        elif which_neighbors == 1:
            return self.get_var_or_file(start_neighbors_from_pdb, params)
        else:
            raise Exception

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

        print >> sys.stderr, 'start to retrieve dists'
        dists = self.get_var_or_file(general_distance, params, recalculate, True, False, False)
        print >> sys.stderr, 'finish to retrieve dists'
        length = len(seq)
        avg_deg = self.get_param(params, 'avg_deg')
        edges = []

        print >> sys.stderr, 'start neighbor'

        effective_length = None

        import wrapper
        msa = self.get_var_or_file(wrapper.my_msa_obj_wrapper, params, False, False, False)

        assert length == msa.get_alignment_length()
        if self.get_param(params, 'which_dist') == 3:
            effective_length = 0
            r = self.get_param(params, 'psicov_gap')
            for i in range(length):
                if float(msa.get_column(i).count('-')) / len(msa) < r:
                    effective_length += 1
        else:
            effective_length = length


        rank_cutoff = ((effective_length * 1.0) * avg_deg)

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

        print >> sys.stderr, 'end neighbor'

        
        return edges


class protease_stone_mutation_list_all(wrapper.obj_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        params.set_param('uniprot_id', global_stuff.PROTEASE_GENE_NAME_TO_USE)


        short_seq = SeqIO.read(global_stuff.PROTEASE_SHORT_FILE, 'fasta')
        long_seq = self.get_var_or_file(dW, params)
        short_to_long, long_to_short = helper.get_alignment_mapping(short_seq.seq.tostring(), long_seq.seq.tostring())
        

        mutation_list = []
        f = open(global_stuff.PROTEASE_MUT_FILE, 'r')
        f.next()
        for line in f:
            s = line.strip().split('\t')
            wild = s[0]
            short_pos = int(s[1]) - 1
            mut = s[2]
            name = global_stuff.PROTEASE_GENE_NAME_TO_USE
            cls = s[6].strip()
            try:
                long_pos = short_to_long[short_pos]
                assert long_seq[long_pos] == wild
                mutation = [name, long_pos, wild, mut, cls]
                mutation_list.append(mutation)
            except Exception, err:
                print err
                pdb.set_trace()
                pass

        return mutation_list

class filtered_mutation_list(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        which_dataset = params.get_param('which_dataset')
        if which_dataset == 'cosmic':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | current_cosmic_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'their_cosmic':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | their_cosmic_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'humvar':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | humvar_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'saapdb':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | saapdb_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'p53':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | p53_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'p53_stone':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | p53_stone_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'protease_stone':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | protease_stone_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'lys_stone':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | lys_stone_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'hemo_stone':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | hemo_stone_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'reverse_stone':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | reverse_stone_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'humsavar':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | humsavar_mutation_list_all.get_all_keys(params, self)
        elif which_dataset == 'cbs':
            return set(['which_dataset', 'protein_list_file', 'mut_freq']) | cbs_mutation_list_all.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return False

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        which_dataset = self.get_param(params, 'which_dataset')

        if which_dataset == 'cosmic':
            mutation_list = self.get_var_or_file(current_cosmic_mutation_list_all, params)
        elif which_dataset == 'humvar':
            mutation_list = self.get_var_or_file(humvar_mutation_list_all, params)
        elif which_dataset == 'their_cosmic':
            mutation_list = self.get_var_or_file(their_cosmic_mutation_list_all, params)
        elif which_dataset == 'saapdb':
            mutation_list = self.get_var_or_file(saapdb_mutation_list_all, params)
        elif which_dataset == 'p53':
            mutation_list = self.get_var_or_file(p53_mutation_list_all, params)
        elif which_dataset == 'p53_stone':
            mutation_list = self.get_var_or_file(p53_stone_mutation_list_all, params)
        elif which_dataset == 'protease_stone':
            mutation_list = self.get_var_or_file(protease_stone_mutation_list_all, params)
        elif which_dataset == 'reverse_stone':
            mutation_list = self.get_var_or_file(reverse_stone_mutation_list_all, params)
        elif which_dataset == 'lys_stone':
            mutation_list = self.get_var_or_file(lys_stone_mutation_list_all, params)
        elif which_dataset == 'hemo_stone':
            mutation_list = self.get_var_or_file(hemo_stone_mutation_list_all, params)
        elif which_dataset == 'humsavar':
            mutation_list = self.get_var_or_file(humsavar_mutation_list_all, params)
        elif which_dataset == 'cbs':
            mutation_list = self.get_var_or_file(cbs_mutation_list_all, params)

        #return mutation_list

        num = 0
        num_add = 0

        filtered_mutations = []



        #m = self.get_var_or_file(general_start_to_pdb, params)

        num_no_neighbor = 0
        mut_freq = params.get_param('mut_freq')
        protein_list = helper.get_file_string_set(self.get_param(params, 'protein_list_file'))

        #g = self.get_var_or_file(general_start_to_pdb_definitive, params)

        for mutation in mutation_list:

            protein_name = mutation[0]

            if num % 50 == 0:
                print num_add, num
            num += 1

            pos = mutation[1]
            wild_res = mutation[2]
            mut_res = mutation[3]
            #change = sum(mutation[4:]) / 8.0
            
            try:

                self.set_param(params, 'uniprot_id', protein_name)
                seq = self.get_var_or_file(dW, params)
                #filtered_mutations.append(mutation)
                
                assert seq[pos] == wild_res


                msa = self.get_var_or_file(wrapper.my_msa_obj_wrapper, params)
                
                col = msa.get_column(pos)
                
                
                n = self.get_var_or_file(general_neighbors_w_weight_w, params)
                deg = len(n[pos])

                if deg == 0:
                    num_no_neighbor += 1
                if col.count(mut_res) > mut_freq and deg > 0:
                
                    filtered_mutations.append(mutation)
                    num_add += 1
                    print num_add, len(mutation_list), len(msa)
                    #import wrapper
                    #self.set_param(params, 'uniprot_id', protein_name)
                    #seq = self.get_var_or_file(dW, params)
                    #msa = self.get_var_or_file(wrapper.my_msa_obj_wrapper, params)
                    #n = self.get_var_or_file(general_neighbors_w_weight_w, params)
                    #deg = len(n[pos])

                    #assert seq[pos] == wild_res
                    #if protein_name in g and seq[pos] == wild_res and msa.get_column(pos).count(mut_res) > 4:
                    #if True:
                    #    filtered_mutations.append(mutation)
                    #    num_add += 1
            except:
                    pass
        print 'num_add: ', num_add, num_no_neighbor

        return filtered_mutations


class cbs_mutation_list_all(wrapper.obj_wrapper):

    def whether_to_override(self, object_key):
        return True

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()


    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):    
        f = open(global_stuff.CBS_MUTATIONS, 'r')
        mutation_list = []
        for line in f:
            s = line.strip().split(',')
            wild = s[0]
            pos = int(s[1])-1
            mut = s[2]
            try:
                one_level = float(s[3])
            except:
                one_level = 0.0
            try:
                two_level = float(s[4])
            except:
                two_level = 0.0
            mutation = ['CBS', pos, wild, mut, one_level, two_level]
            mutation_list.append(mutation)
        return mutation_list


class hemo_stone_mutation_list_all(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):    
        f = open(global_stuff.HEMO_MUT_FILE, 'r')
        f.next()
        mutation_list = []
        params.set_param('uniprot_id', global_stuff.HEMO_GENE_NAME_TO_USE)
        seq = self.get_var_or_file(dW, params)
        for line in f:
            s = line.strip().split('\t')
            wild = s[0]
            mut = s[2]
            pos = int(s[1])-1
            cls = s[6].strip()
            assert seq[pos] == wild
            mutation = [global_stuff.HEMO_GENE_NAME_TO_USE, pos, wild, mut, cls]
            mutation_list.append(mutation)
        return mutation_list

class p53_stone_mutation_list_all(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(global_stuff.STONE_P53_FILE)
        mutation_list = []
        f.next()
        for line in f:
            s = line.strip().split('\t')
            name = 'TP53'
            pos = int(s[1]) - 1
            wild = s[0]
            mut = s[2]
            freq = int(s[6])
            mutation = [name, pos, wild, mut, freq]
            mutation_list.append(mutation)
        return mutation_list


# one wrapper for each raw dataset.  no parameters needed
class their_cosmic_mutation_list_all(wrapper.obj_wrapper):

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        f = open(global_stuff.their_mutation_file, 'r')
        mutations = []
        total = 0
        bad = 0
        for line in f:
            try:
                s = line.strip().split('\t')
                name = s[1]
                pos_raw = s[6]
                wild_char = pos_raw[0]
                mut_char = pos_raw[-1]
                pos = int(pos_raw[1:-1]) - 1
                mut_count = int(s[7])
                site_count = int(s[8])
                gene_count = int(s[10])
                func = int(s[15] == '1')
                protein_binding = int(s[16] == '1')
                dna_binding = int(s[17] == '1')
                mol_binding = int(s[18] == '1')
                specificity_score = float(s[14])
                conservation_score = float(s[13])
                their_score = float(s[12])
                params.set_param('uniprot_id', name)
                #import wc
                #seq = wc.get_stuff(dW, params)
                #assert seq[pos] == wild_char
                mutation = [name, pos, wild_char, mut_char, mut_count, site_count, gene_count, func, protein_binding, dna_binding, mol_binding, their_score, conservation_score, specificity_score]
                mutations.append(mutation)
            except Exception, err:
                bad += 1
                print name
                print err
            total += 1
        print 'bad: ', bad, ' total: ', total

        return mutations


    

def get_protein_info(protein_list, info_file, params):
    import wc
    f = open(protein_list, 'r')
    g = open(info_file, 'w')
    i = 0
    for line in f:
        name = line.strip()
        print >> sys.stderr, name,i
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
            
            def get_blastp_msa_length():
                params.set_param('which_blast', 1)
                msa = wc.get_stuff(my_msa_obj_wrapper, params, False, False, False)
                return [str(len(msa))]

            def get_delta_blast_msa_length():
                params.set_param('which_blast', 2)
                assert params.get_param('which_msa') == 2
                import wrapper
                msa = wc.get_stuff(wrapper.my_msa_obj_wrapper, params, False, False, False)
                return [str(len(msa))]

            def get_num_mutations():
                mutations = wc.get_stuff(protein_mutation_list, params, False, False, False)
                return [str(len(mutations))]

            info = []
            which_info = [get_name, get_seq_length, get_delta_blast_msa_length, get_num_mutations]
            #which_info = [get_name, get_seq_length, get_msa_length]
            for which in which_info:
                info = info + which()
            
            g.write(string.join(info,sep=',') + '\n')
            print global_stuff.time_total
        except:

            print >> sys.stderr, 'fail'
    f.close()
    g.close()
    print global_stuff.time_total


def get_every_site_info(params, protein_list, dist_file):
    import wc
    f = open(dist_file,'w')

    protein_list = protein_list
    h = open(protein_list)

    k=0

    for line in h:
        protein_name = line.strip()
        print >> sys.stderr, protein_name, k
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

        print >> sys.stderr, helper.get_overlap(n1,n2)

    


def get_mutation_info(out_file, params):
    import wc

    
    l = wc.get_stuff(filtered_mutation_list, params, False, False, False, False)

    f = open(out_file, 'w')

    avg_degs = [1,2,3,4]
    
    i = 0
    pdb.set_trace()
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
            seq = wc.get_stuff(dW, params)
            print seq[mutation[1]], mutation[2], mutation
            assert seq[mutation[1]] == mutation[2]
            
            msa = wc.get_stuff(wrapper.my_msa_obj_wrapper, params, False, False, False)
            col = msa.get_column(mutation[1])
            return [str(col.count(mutation[2]))]

        def get_avg_change():
            return [str(sum(mutation[4:]) / 8.0)]

        def get_wild_aa():
            seq = wc.get_stuff(dW, params)
            pos = mutation[1]
            wild = mutation[2]
            #try:
            #    assert seq[pos] == wild
            #except Exception, err:
            #    pdb.set_trace()
            #    print err
            return [wild]

        def get_mut_aa():
            seq = wc.get_stuff(dW, params)
            pos = mutation[1]
            mut = mutation[3]

            return [mut]


        def get_wild_category_num():
            msa = wc.get_stuff(wrapper.my_msa_obj_wrapper, params, False, False, False)
            col = msa.get_column(mutation[1])
            wild_res_cat = global_stuff.aa_to_class[mutation[2]]
            count = 0
            for aa in col:
                if global_stuff.aa_to_class[aa] == wild_res_cat:
                    count += 1
            return [str(count)]

        def get_mut_category_num():
            msa = wc.get_stuff(wrapper.my_msa_obj_wrapper, params, False, False, False)
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
            msa = wc.get_stuff(wrapper.my_msa_obj_wrapper, params, False, False, False)
            col = msa.get_column(mutation[1])
            return [str(col.count(mutation[3]))]

        def get_msa_length():
            msa = wc.get_stuff(wrapper.my_msa_obj_wrapper, params, False, False, False)
            return [str(len(msa))]

        def get_cosmic_info():
            return [str(mutation[4]), str(mutation[5]), str(mutation[6])]

        def get_whether_bad():
            return [str(mutation[4])]

        def get_deg_pdb():
            cos = [5.0,6.0,7.0,8.0]
            pos = mutation[1]
            ans = []
            for co in cos:
                params.set_param('co', co)
                neighbors = wc.get_stuff(general_neighbors_w_weight_w, params)
                deg = len(neighbors[pos])
                ans.append(str(deg))
            return ans

        def get_p53_stone_count():
            return [str(mutation[4])]

        info = []
        if i % 50 == 0:
            print >> sys.stderr, get_name(), i

        which_dataset = params.get_param('which_dataset')
        #if which_dataset == 'cosmic':

        #    which_info = [get_name, get_pos, get_wild_aa, get_mut_aa, get_cosmic_info]

        #elif which_dataset == 'humvar':
        #    which_info = [get_name, get_pos, get_wild_num, get_mut_num, get_deg, get_whether_bad, get_msa_length]

        #elif which_dataset == 'their_cosmic':
        print i
        if i == 2303:
            pdb.set_trace()
        i += 1
        #which_info = [get_name, get_pos, get_wild_num, get_mut_num, get_msa_length, get_wild_num, get_mut_num, get_deg_pdb]
        which_info = [get_name, get_pos, get_wild_num, get_mut_num, get_msa_length, get_deg_pdb]
        try:
            for which in which_info:
                info = info + which()

            f.write(string.join(info,sep=',') + '\n')
        except Exception, err:
            pdb.set_trace()
            print >> sys.stderr, err
            pass
        

    f.close()

# outputs file.  could be roc file input, or other input.  one argument is function that assigns a number to each mutation.  the function determines the output
def get_output_obj(params, l, use_neighbor, ignore_pos, max_neighbors, num_trials, pseudo_total, sim_f, norm_f, mut_to_num_f, to_neighbor_p_value):
    import wc
    #params.set_param('protein_list_file', protein_list_file)
    #l = wc.get_stuff(filtered_mutation_list_given_protein_list, params)
    i = 0
    scores = []
    labels = []
    bad_count = 0



    for mutation in l:

        if i%1 == 0:
            try:
                from mpi4py import MPI

                comm = MPI.COMM_WORLD
                rank = comm.Get_rank()
            except:
                rank = 0
            if i % 1 == 0:
                print i, 'calculating, rank: ', rank


        i += 1

        try:



            score = helper.predict_position_energy_weighted(params, mutation, use_neighbor, ignore_pos, max_neighbors, num_trials, pseudo_total, sim_f, to_neighbor_p_value)
            scores.append(score)

            labels.append(mut_to_num_f(mutation) + mutation)
            
        except Exception, err:
            bad_count += 1
        
            import traceback
            for frame in traceback.extract_tb(sys.exc_info()[2]):
                fname,lineno,fn,text = frame
                print "Error in %s on line %d" % (fname, lineno)
            print >> sys.stderr, err, mutation, bad_count




    assert len(scores) == len(labels)

    normed_scores = norm_f(scores)
    ans = [ None for i in range(len(scores))]
    for i in range(len(scores)):
        ans[i] = [normed_scores[i]] + labels[i]


    return ans

    #helper.write_mat(ans, out_file)


    
    

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
            print >> sys.stderr, '      ', i


        i += 1
        try:
            score = helper.predict_position_energy_weighted(params, recalculate, mutation, use_neighbor, ignore_pos, max_neighbors, weighted, num_trials, pseudo_total, sim_f)
            ans.append([score, mutation[4]])
        except Exception, err:
            bad_count += 1
            print >> sys.stderr, 'failed to get score', err, bad_count

    helper.write_mat(ans, out_file)



from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.PDB import Polypeptide
from Bio.Blast import NCBIXML
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO



class fW(wrapper.file_wrapper):

    def whether_to_override(self, object_key):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['pdb'])

    def get_folder(self, object_key):
        return global_stuff.BIN_FOLDER + 'fWs' + '/'

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        pdb_file_name = self.get_param(params, 'pdb')
        pdbl = Bio.PDB.PDBList()
        pdbl.retrieve_pdb_file(pdb_file_name, pdir=self.get_holding_folder())
        subprocess.call(['mv', self.get_holding_folder() + string.lower('pdb'+pdb_file_name+'.ent'), self.get_holding_location()])

        return open(self.get_holding_location(), 'r')


class cW(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):
    """
    returns chain of pdb chain, which is a list of residue objects
    this is the object i find neighbors with
    """

    def whether_to_override(self, object_key):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['chain']) | fW.get_all_keys(params, self)


    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        f = self.get_var_or_file(fW, params, recalculate, to_pickle)
        structure = Bio.PDB.PDBParser().get_structure(self.get_param(params, 'pdb'), f)
        chain = Bio.PDB.PPBuilder().build_peptides(structure[0][self.get_param(params, 'chain')])

        to_return = []
        for chain_frag in chain:
            to_return = to_return + chain_frag


        real_to_return = []

        for i in range(len(to_return)):

            try:
                Polypeptide.three_to_one(to_return[i].resname)
            except Exception, err:
                print err
                print to_return[i], i
                #pdb.set_trace()
                print 'weird AA', to_return[i].resname
            else:
                real_to_return.append(to_return[i])

        return real_to_return

class pdb_chain_seq(wrapper.obj_wrapper, wrapper.by_pdb_folder_wrapper):

    def whether_to_override(self, object_key):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        return cW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        chain = self.get_var_or_file(cW, params)
        pdb_name = self.get_param(params, 'pdb')
        chain_name = self.get_param(params, 'chain')
        raw_string = [Polypeptide.three_to_one(res.resname) for res in chain]
        
        raw_string = string.join(raw_string, sep = '')
        ans = SeqRecord(Seq(raw_string), id = pdb_name + '_' + chain_name)

        return ans

class gene_name_to_refseq(wrapper.obj_wrapper):


    def whether_to_override(self, object_key):
        return False

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(global_stuff.their_mutation_file, 'r')
        m = {}
        for line in f:
            s = line.strip().split('\t')
            gene_name = s[1].strip()
            refseq = s[4].strip()
            m[gene_name] = refseq

        return m
    

class refseq_to_pdb(wrapper.obj_wrapper):
    """
    maps uniprot id to pdb, chain tuple
    """

    def whether_to_override(self, object_key):
        return False
    
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        r_to_u = self.get_var_or_file(refseq_to_uniprot, params)
        u_to_p = self.get_var_or_file(uniprot_to_pdb, params)
        m = {}
        i = 0
        for refseq in r_to_u:
            print i
            i += 1
            for uniprot in r_to_u[refseq]:
                try:
                    for pc in u_to_p[uniprot]:
                        try:
                            m[refseq].append(pc)
                        except:
                            m[refseq] = [pc]
                except:
                    pass
        return m
            

class refseq_to_uniprot(wrapper.obj_wrapper):
    """
    assuming that there is only one refseq file ever, so that i don't have to store the specific refseq file as a key
    """

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(global_stuff.refseq_to_uniprot_file, 'r')
        m = {}
        for line in f:
            s = line.strip().split(',')
            refseq = s[0]
            uniprot = s[1]
            try:
                m[refseq].append(uniprot)
            except:
                m[refseq] = [uniprot]
        return m

# this can be used by different datasets
class uniprot_to_pdb(wrapper.obj_wrapper):
    """
    one to many map from uniprot to pdb.  
    """
    def whether_to_override(self, object_key):
        return False
    
    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(global_stuff.pdb_to_uniprot_file, 'r')
        m = {}
        for line in f:
            s = line.strip().split(',')
            pdb_name = s[0].upper()
            chain_letter = s[1].upper()
            uniprot_id = s[2]
            to_add = (pdb_name, chain_letter)
            try:
                m[uniprot_id].append(to_add)
            except:
                m[uniprot_id] = [to_add]
        return m

class refseq_to_pdb_definitive(wrapper.obj_wrapper):
    """
    contains all the refseqs which were in the cosmic file from sander.  but when looking at his file, some of corresponding gene names i may not have the file for
    """

    def whether_to_override(self, object_key):
        return False

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        m = self.get_var_or_file(refseq_to_pdb, params)
        definitive = {}
        i = 1
        for refseq in m:
            print i
            i += 1
            longest = -1
            for pc in m[refseq]:
                pdb_name = pc[0]
                chain_letter = pc[1]
                self.set_param(params, 'pdb', pdb_name)
                self.set_param(params, 'chain', chain_letter)
                chain = self.get_var_or_file(pdb_chain_seq, params)
                if len(chain) > longest:
                    definitive[refseq] = pc
                    longest = len(chain)
        return definitive

class saapdb_to_pdb(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return uniprot_to_pdb.get_all_keys(params, self)

    def whether_to_override(self, object_key):
        return False

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        saapdb_genes = helper.get_file_string_set(global_stuff.saapdb_genes_file)
        u_to_p = self.get_var_or_file(uniprot_to_pdb, params)
        m = {}
        for u in u_to_p:
            if u in saapdb_genes:
                m[u] = u_to_p[u]
        return m

# assume there is mapping from start id to possible list of pdb's.  this can then be used by general_start_to_pdb (code already there is badly organized special case
# keys contain only id's for which there are pdb's
class general_start_to_pdb(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):

        which_dataset = params.get_param('which_dataset')
        if which_dataset == 'saapdb':
            return set(['which_dataset']) | saapdb_to_pdb.get_all_keys(params, self)
        elif which_dataset == 'their_cosmic':
            pass
            #return set(['which_dataset']) | saapdb_to_pdb.get_all_keys(params, self)
        elif which_dataset == 'humsavar':
            return humsavar_start_to_pdb.get_all_keys(params, self)


    def whether_to_override(self, object_key):
        return False

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        which_dataset = self.get_param(params, 'which_dataset')
        if which_dataset == 'saapdb':
            return self.get_var_or_file(saapdb_to_pdb, params)
        elif which_dataset == 'humsavar':
            return self.get_var_or_file(humsavar_start_to_pdb, params)





class humsavar_mutation_list_all(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(global_stuff.HUMSAVAR_FILE, 'r')
        mutation_list = []
        for line in f:
            s = line.strip().split()
            start = s[0]
            uniprot = s[2]
            pos = int(s[4])
            wild = s[5]
            mut = s[7]
            label = s[8]
            mutation = [start, pos, wild, mut, label]
            mutation_list.append(mutation)
        return mutation_list


class humsavar_start_to_pdb(wrapper.obj_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set()

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        f = open(global_stuff.HUMSAVAR_FILE, 'r')
        m = {}
        u_to_p = self.get_var_or_file(uniprot_to_pdb, params)
        for line in f:
            s = line.strip().split()
            start = s[0]
            uniprot = s[2]
            pos = int(s[4])
            wild = s[5]
            mut = s[7]
            label = s[8]

            if uniprot in u_to_p:
                m[start] = u_to_p[uniprot]
            
        return m





class general_start_to_pdb_definitive(wrapper.obj_wrapper):
    """
    only contains mapping of gene names for which i have the sequence
    """

    @classmethod
    def get_all_keys(cls, params, self=None):
        which_data = params.get_param('which_dataset')
        if which_data == 'p53':
            return set(['which_dataset'])
        elif which_data == 'cbs':
            return set(['which_dataset'])
        elif which_data == 'p53_stone':
            return set(['which_dataset'])
        elif which_data == 'protease_stone':
            return set(['which_dataset'])
        elif which_data == 'lys_stone':
            return set(['which_dataset'])
        elif which_data == 'reverse_stone':
            return set(['which_dataset'])
        elif which_data == 'hemo_stone':
            return set(['which_dataset'])
        elif which_data == 'their_cosmic':
            return set(['which_dataset']) | gene_name_to_refseq.get_all_keys(params, self) | refseq_to_pdb_definitive.get_all_keys(params, self)
        else:
            return (set(['which_dataset']) | general_start_to_pdb.get_all_keys(params, self) | pdb_chain_seq.get_all_keys(params, self)) - set(['pdb', 'chain'])

    def whether_to_override(self, object_key):
        return True

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        which_dataset = self.get_param(params, 'which_dataset')
        # this is the special case.  otherweise, do all of the computing here
        if which_dataset == 'p53':
            m = {'TP53':('1TUP','A')}
            return m
        elif which_dataset == 'p53_stone':
            m = {'TP53':('1TUP','A')}
            return m
        elif which_dataset == 'protease_stone':
            m = {global_stuff.PROTEASE_GENE_NAME_TO_USE:('5HVP','A')}
            return m
        elif which_dataset == 'reverse_stone':
            m = {global_stuff.REVERSE_GENE_NAME_TO_USE:('3KJV','B')}
            return m
        elif which_dataset == 'lys_stone':
            m = {'LYSOZYME':('102L','A')}
            return m
        elif which_dataset == 'hemo_stone':
            m = {global_stuff.HEMO_GENE_NAME_TO_USE:('3EOK','B')}
            return m
        elif which_dataset == 'cbs':
            m = {'CBS':('1JBQ','A')}
            return m
        elif which_dataset == 'their_cosmic':
            # only add gene name for which i have the sequence
            theirs_i_have = helper.get_file_string_set(global_stuff.cosmic_genes_file)


            
            gn_to_r = self.get_var_or_file(gene_name_to_refseq, params)
            r_to_p = self.get_var_or_file(refseq_to_pdb_definitive, params)
            m = {}
            for gn in gn_to_r:
                if gn_to_r[gn] in r_to_p and gn in theirs_i_have:
                  m[gn] = r_to_p[gn_to_r[gn]]
            return m
        else:
            m = self.get_var_or_file(general_start_to_pdb, params)
            definitive = {}
            i = 1
            for s in m:
                print i
                i += 1
                longest = -1
                for pc in m[s]:
                    pdb_name = pc[0]
                    chain_letter = pc[1]
                    self.set_param(params, 'pdb', pdb_name)
                    self.set_param(params, 'chain', chain_letter)
                    try:
                        chain = self.get_var_or_file(pdb_chain_seq, params)
                        if len(chain) > longest:
                            definitive[s] = pc
                            longest = len(chain)
                    except:
                        pass
            return definitive            
    
class query_to_pdb_chain_maps(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):
    """
    tuple of uniprot to pdb and pdb to uniprot position map
    """

    @classmethod
    def get_all_keys(cls, params, self=None):
        return ( gene_name_to_refseq.get_all_keys(params, self) | uniprot_to_pdb.get_all_keys(params, self) | pdb_chain_seq.get_all_keys(params, self) | dW.get_all_keys(params, self) ) \
               - set(['pdb', 'chain'])

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        uniprot_id = self.get_param(params, 'uniprot_id')

        m = self.get_var_or_file(general_start_to_pdb_definitive, params)

        pdb_name, chain_letter = m[uniprot_id]
        self.set_param(params, 'pdb', pdb_name)
        self.set_param(params, 'chain', chain_letter)
        pdb_chain = self.get_var_or_file(pdb_chain_seq, params)
        uniprot_chain = self.get_var_or_file(dW, params)
        uniprot_seq = uniprot_chain.seq.tostring()
        pdb_seq = pdb_chain.seq.tostring()

        pdb.set_trace()

        return helper.get_alignment_mapping(uniprot_seq, pdb_seq)


        
        matrix = matlist.blosum62
        gap_open = -10
        gap_extend = -0.5

        alns = pairwise2.align.globalds(uniprot_seq, pdb_seq, matrix, gap_open, gap_extend)
        top_aln = alns[0]
        aln_uniprot, aln_pdb, score, begin, end = top_aln
        uniprot_i = 0
        pdb_i = 0
        aln_length = len(aln_uniprot)

        uniprot_to_pdb_map = {}
        pdb_to_uniprot_map = {}



        for i in range(aln_length):
            if aln_pdb[i] != '-' and aln_uniprot[i] != '-':
                uniprot_to_pdb_map[uniprot_i] = pdb_i
                pdb_to_uniprot_map[pdb_i] = uniprot_i
                
            if aln_pdb[i] != '-':
                pdb_i += 1
            if aln_uniprot[i] != '-':
                uniprot_i += 1

        return uniprot_to_pdb_map, pdb_to_uniprot_map


class chain_distances(wrapper.mat_obj_wrapper, wrapper.by_pdb_folder_wrapper):

    def whether_to_override(self, object_key):
        return False

    @classmethod
    def get_all_keys(cls, params, self=None):
        return cW.get_all_keys(params, self)

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        print 'CALCULATING DISTANCES!!!!', params.get_param('pdb'), params.get_param('chain')
        residues = self.get_var_or_file(cW, params, recalculate, True)
        rep_atoms = [helper.get_representative_atom(res) for res in residues]
        num_res = len(residues)
        dists = [[-1 for i in range(num_res)] for j in range(num_res)]
        for i in range(num_res):
            for j in range(num_res):
                try:
                    dists[i][j] = rep_atoms[i] - rep_atoms[j]
                except Exception as e:
                    print 'ERROR: distance fail', self.params, i, j, residues[i].child_dict.keys(), residues[j].child_dict.keys()
                    dists[i][j] = -1
        return dists

# if is a by_uniprot_id_wrapper, that means gene_name/uniprot_id are specified in params
class start_neighbors_from_pdb(wrapper.int_float_tuple_mat_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @classmethod
    def get_all_keys(cls, params, self=None):
        return set(['co']) | (query_to_pdb_chain_maps.get_all_keys(params, self) | query_to_pdb_chain_maps.get_all_keys(params, self) | general_start_to_pdb_definitive.get_all_keys(params, self)) 

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        # assume uniprot_id (actually gene name) is specified

        uniprot_id = self.get_param(params, 'uniprot_id')
        u_to_p, p_to_u = self.get_var_or_file(query_to_pdb_chain_maps, params)
        pdb_name, chain_letter = self.get_var_or_file(general_start_to_pdb_definitive, params)[uniprot_id]
        self.set_param(params, 'pdb_name', pdb_name)
        self.set_param(params, 'chain_letter', chain_letter)
        dists = self.get_var_or_file(chain_distances, params)
        uniprot_seq = self.get_var_or_file(dW, params)
        pdb_seq = self.get_var_or_file(pdb_chain_seq, params)
        cutoff = self.get_param(params, 'co')
        assert len(dists) == len(pdb_seq)
        neighbors = [ [] for i in range(len(uniprot_seq)) ]
        for i in range(len(uniprot_seq)):

            try:
                pdb_pos = u_to_p[i]
            except:
                pass
            else:
                for j in range(len(pdb_seq)):
                    if dists[pdb_pos][j] != -1 and dists[pdb_pos][j] < cutoff and pdb_pos != j:
                        try:
                            neighbors[i].append((p_to_u[j], dists[pdb_pos][j]))
                        except:
                            pass

        return neighbors
        
    def whether_to_override(self, object_key):
        return True


def write_genes_mutation_list(p, outfile):
    import wc
    mutation_list = wc.get_stuff(filtered_mutation_list, p)
    gene_list = set([mutation[0] for mutation in mutation_list])
    helper.write_iterable_vert(gene_list, outfile)

    
