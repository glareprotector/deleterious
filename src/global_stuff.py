'''
Created on Mar 9, 2012

@author: glareprotector
'''
import string
import csv
import string
import re
import math

import pdb

def get_param():
    import param
    p = param.param({'ev':.05, 'uniprot_id':'test', 'avg_deg':3, 'n_cutoff':0, 'f_cutoff':15, 'which_msa':0, 'which_weight':0, 'which_dist':0, 'pseudo_c':0.1, 'which_blast':1, 'blmax':700})
    return p


cosmic_or_humvar = 'cosmic'

#home = '/home/fw27/d/deleterious/'
home = '/mnt/work/fultonw/deleterious/'
home = '/home/fw27/d/deleterious/'
#home = '/mnt/work/fultonw/deleterious/'
#base_folder = '../data/proteins/humvar_from_orchestra/'
#base_folder = '/mnt/work/fultonw/deleterious/data/proteins/humvar/'
#base_folder = '/mnt/work/fultonw/scratch/'
base_folder = '/mnt/work/fultonw/scratch_cosmic/'
#base_folder = '/home/fw27/d/deleterious/data/proteins/humvar_from_orchestra/'
#base_folder = '/mnt/work/fultonw/deleterious/data/proteins/humvar/'
#base_folder = '/home/fw27/d/deleterious/data/proteins/cosmic/'




BIN_FOLDER = home + 'data/bin/'
HOLDING_FOLDER = home + 'data/holding_folder/'
lock_folder = home + 'lock_folder/'
process_folder = home + 'process_files/'
polyphen_msa_directory = home + 'data/polyphen-2.2.2/precomputed/alignments/'
data_folder = home + 'data/'


# path to files
all_seqs_file = '../data/human-2011_12.seq'
neutral_mutations_file = '../data/humvar-2011_12.neutral.pph.input'
deleterious_mutations_file = '../data/humvar-2011_12.deleterious.pph.input'
cosmic_raw_data_folder = data_folder + 'fasta/'






# path to programs
MUSCLE_PATH = '/home/fw27/d/deleterious/muscle3.8.31_i86linux64'
BLAST_PATH = '/mnt/work/fultonw/deleterious/blast/ncbi-blast-2.2.26+/bin/psiblast'
BLAST_PATH = '/home/fw27/d/deleterious/bin/psiblast'
#BLASTP_PATH = '/mnt/work/fultonw/deleterious/blast/ncbi-blast-2.2.26+/bin/blastp'
BLASTP_PATH = '/home/fw27/d/deleterious/bin/blastp'
#BLASTDB_PATH = '/mnt/work/fultonw/nr/'
#BLASTDB_PATH = 'nr/nr'
BLASTDB_PATH = '/groups/shared_databases/blastdb/nr'
MIP_PATH = 'MIp_wrapper.pl'


# random constants
query_gi_number = '123456789123456789'
proc_id = 0
whether_to_look_at_whether_to_override = False
to_reindex = True
recalculate = False
recalculate_nodewise_loss_f = True
metric_cutoffs = [1,2,3,4,5,6,7,8,9]
aa_to_num = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
aa_to_class = {'A':0, 'G':0, 'V':0, 'I':1, 'L':1, 'P':1, 'F':1, 'Y':2, 'T':2, 'M':2, 'S':2, 'C':2, 'H':3, 'N':3, 'E':3, 'W':3, 'R':4, 'K':4, 'D':5, 'Q':5, 'X':6, 'B':6, 'Z':6, '-':7}
ignore_aas = ['-','X', 'B', 'Z']
q = 20
#aa_to_num = {'A':2,'M':1,'C':0}
#q=3
