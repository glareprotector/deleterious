
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
from Bio import AlignIO

import math
import subprocess
import string
import os
import random
import pdb
import re






class bW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):
        # all this will do is read specified fasta file, and save it in my format
        protein_name = self.get_param(params, 'uniprot_id')
        protein_folder = global_stuff.PROTEIN_BASE_FOLDER + protein_name + '/'
        existing_seq_file = protein_folder + 'seq'

        old_seq = SeqIO.read(open(existing_seq_file,'r'), 'fasta')

        seq_record = Bio.SeqRecord.SeqRecord(old_seq)
        SeqIO.write(old_seq, self.get_holding_location(), 'fasta')
        return open(self.get_holding_location(),'r')

class dW(wrapper.obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle, to_filelize = False, always_recalculate = False, old_obj = None):

        seq_file_handle = self.get_var_or_file(bW, params, recalculate, True, False)
        asdf = SeqIO.read(seq_file_handle, 'fasta')
        return asdf

# blast results file wrapper(xml format)
class adW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):



    def whether_to_override(self, location):
        
        #if the file size is too small, we know there was something wrong
        import os
        if os.path.getsize(location) < 1:
            return True

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        seq_records = []
        f = self.get_var_or_file(bW, params, recalculate, False, False, False)
        query = SeqIO.parse(f, 'fasta')
        seq_records.append(query)

        print 'RUNNING BLAST!!!!!!!'
        psi_blast_cline = NcbipsiblastCommandline(cmd = global_stuff.BLAST_PATH, outfmt = 5, query = '\''+f.name+'\'', db = global_stuff.BLASTDB_PATH, out = self.get_holding_location())
        #pdb.set_trace()
        subprocess.call(str(psi_blast_cline), shell=True, executable='/bin/bash')

        return open(self.get_holding_location())


# processsed blast results format(bunch of sequences in fasta in 1 file)
class aeW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        # parse blast xml file, then do processing

        blast_xml_handle = self.get_var_or_file(adW, params, recalculate, False, False, False)
        record = NCBIXML.read(blast_xml_handle)
        seen = set()
        seq_records = []
        # add the query sequence, and have a set so that only add each sequence once
        query = self.get_var_or_file(dW, params, recalculate, True, False, False)
        query.id = 'QUERY'
        seq_records.append(query)
        seen.add(query.seq.tostring())
        # add high scoring pairs in alignments with sufficiently low evalue that have not been seen
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < self.get_param(params, 'ev') and not hsp.sbjct in seen:
                    seq_records.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(hsp.sbjct), id = alignment.hit_id))
        # write hits to fasta file
        output_handle = open(self.get_holding_location(), 'w')
        SeqIO.write(seq_records, output_handle, 'fasta')
        print 'WROTE ', self.get_holding_location()
        return output_handle


# gets the result of muscle
class afW(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        msa_input_handle = self.get_var_or_file(aeW, params, recalculate, to_pickle, to_filelize, always_recalculate)
        cline = MuscleCommandline(cmd = global_stuff.MUSCLE_PATH, input = '\''+msa_input_handle.name+'\'', out = self.get_holding_location(), clw = False, maxiters = 2)
        #cline()

        subprocess.call(str(cline), shell=True, executable='/bin/bash')
        return open(self.get_holding_location())


# processed msa output(columns with skips removed)
class agW(wrapper.msa_obj_wrapper, wrapper.by_uniprot_id_wrapper):

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

class easy_to_read_msa_format(wrapper.file_wrapper, wrapper.by_uniprot_id_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):

        f = open(self.get_holding_location(), 'w')
        recalculate = False
        pdb.set_trace()
        msa = self.get_var_or_file(agW, params, recalculate, True, True, False)
        f.write(str(len(msa)) + '\t' + str(msa.get_alignment_length()) + '\n')
        for i in range(len(msa)):
            to_write = ''
            for j in range(msa.get_alignment_length()):
                to_write = to_write + msa[i,j]
            print len(to_write), to_write
            f.write(to_write + '\n')
        return f

class pairwise_dist(wrapper.mat_obj_wrapper, wrapper.by_uniprot_id_wrapper):

    @dec
    def constructor(self, params, recalculate, to_pickle = False, to_filelize = False, always_recalculate = False, old_obj = None):
        msa = self.get_var_or_file(agW, params, recalculate, False, False, False)
        dists = [ [0 for i in range(msa.get_alignment_length())] for j in range(msa.get_alignment_length())]
        for i in range(msa.get_alignment_length()):
            for j in range(msa.get_alignment_length()):
                x_no_skip = ''
                y_no_skip = ''
                for k in range(len(msa)):
                    if msa[k][i] != '-' and msa[k][j] != '-':
                        x_no_skip += msa[k][i]
                        y_no_skip += msa[k][j]
                dists[i][j] = helper.get_KL(x_no_skip, y_no_skip)
                dists[j][i] = dists[i][j]
        return dists
