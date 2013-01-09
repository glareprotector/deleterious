'''
Created on Mar 9, 2012

@author: glareprotector
'''
import string
import csv
import string
import re
import math
import datetime
import pdb

time_total = datetime.timedelta(0)

cosmic_or_humvar = 'humvar'
orchestra_or_no = 'orchestra'

if orchestra_or_no  == 'orchestra':

    real_home = '/home/fw27/d/deleterious/'
    temp_home = '/tmp/fw27/'
    
    if cosmic_or_humvar == 'humvar':
        real_base_folder = '/home/fw27/d/deleterious/data/proteins/humvar/'
        temp_base_folder = temp_home + 'humvar/'
        remote_base_folder = '/mnt/work/fultonw/scratch/'
    elif cosmic_or_humvar == 'cosmic':
        real_base_folder = '/home/fw27/d/deleterious/data/proteins/cosmic/'
        temp_base_folder = temp_home + 'cosmic/'
        remote_base_folder = '/mnt/work/fultonw/scratch_cosmic/'

    MUSCLE_PATH = '/home/fw27/d/deleterious/muscle3.8.31_i86linux64'
    BLAST_PATH = '/home/fw27/d/deleterious/bin/psiblast'
    BLASTP_PATH = '/home/fw27/d/deleterious/bin/blastp'
    BLASTDB_PATH = '/groups/shared_databases/blastdb/nr'
    DELTABLAST_PATH = '/home/fw27/d/deleterious/src/deltablast'
    CDD_PATH = '/home/fw27/d/deleterious/bin/cdd/cdd_delta'
    MIP_PATH = 'MIp_wrapper.pl'
    PSICOV_PATH = 'psicov'
    HHBLITS_PATH = 'hhblits'
    HHBLITS_DB_PATH = '/home/fw27/d/deleterious/hh/hhdb/nr20_12Aug11'
    HHBLITS_CONVERT_A3M_TO_FASTA = '/home/fw27/d/deleterious/hh/hhsuite-2.0.15-linux-x86_64/lib/hh/scripts/reformat.pl'

elif orchestra_or_no == 'no':

    real_home = '/mnt/work/fultonw/deleterious/'
    temp_home = None

    if cosmic_or_humvar == 'humvar':
        real_base_folder = '/mnt/work/fultonw/scratch/'
        temp_base_folder = None
    elif cosmic_or_humvar == 'cosmic':
        real_base_folder = '/mnt/work/fultonw/scratch_cosmic/'
        temp_base_folder = None
        
    MUSCLE_PATH = '/mnt/work/fultonw/deleterious/muscle/muscle3.8.31_i86linux64'
    BLAST_PATH = '/mnt/work/fultonw/deleterious/blast/ncbi-blast-2.2.26+/bin/psiblast'
    BLASTP_PATH = '/mnt/work/fultonw/deleterious/blast/ncbi-blast-2.2.26+/bin/blastp'
    BLASTDB_PATH = 'nr/nr'

    DELTABLAST_PATH = None
    CDD_PATH = None
    MIP_PATH = 'MIp_wrapper.pl'
    PSICOV_PATH = 'psicov'
    HHBLITS_PATH = 'hhblits'
    HHBLITS_DB_PATH = '/mnt/work/fultonw/deleterious/hh/hhdb/nr20_12Aug11'
    HHBLITS_CONVERT_A3M_TO_FASTA = '/mnt/work/fultonw/deleterious/hh/hhsuite-2.0.15-linux-x86_64/lib/hh/scripts/reformat.pl'



base_folder = real_base_folder
home = real_home

def get_param():
    import param



    p = param.param({'ev':.05, 'uniprot_id':'P16455', 'avg_deg':1, 'n_cutoff':0, 'f_cutoff':15, 'which_msa':2, 'which_weight':1, 'which_dist':3, 'pseudo_c':0.1, 'which_blast':2, 'blmax':999999, 'which_impute':0, 'filter_co':0.35, 'psicov_sep':6, 'psicov_gap':0.5, 'psicov_r':.001, 'psiblast_iter':5, 'hhblits_iter':2})




    return p



def get_home():
    global home
    return home

def get_holding_folder():
    return get_home() + 'data/holding_folder/'

BIN_FOLDER = real_home + 'data/bin/'
HOLDING_FOLDER = get_home() + 'data/holding_folder/'
lock_folder = real_home + 'lock_folder/'
process_folder = real_home + 'process_files/'
log_folder = real_home + 'dump/'
polyphen_msa_directory = real_home + 'data/polyphen-2.2.2/precomputed/alignments/'
data_folder = real_home + 'data/'


# path to files
all_seqs_file = '../data/human-2011_12.seq'
neutral_mutations_file = '../data/humvar-2011_12.neutral.pph.input'
deleterious_mutations_file = '../data/humvar-2011_12.deleterious.pph.input'
cosmic_raw_data_folder = data_folder + 'fasta/'





# random constants
query_gi_number = '123456789123456789'
proc_id = 0
whether_to_look_at_whether_to_override = True
to_reindex = True
recalculate = False
recalculate_nodewise_loss_f = True
metric_cutoffs = [1,2,3,4,5,6,7,8,9]
aa_to_aa = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
aa_to_class = {'A':0, 'G':0, 'V':0, 'I':1, 'L':1, 'P':1, 'F':1, 'Y':2, 'T':2, 'M':2, 'S':2, 'C':2, 'H':3, 'N':3, 'E':3, 'W':3, 'R':4, 'K':4, 'D':5, 'Q':5}
ignore_aas = ['-','X', 'B', 'Z']
aa_to_num = aa_to_class
q = max(aa_to_num.values())+1
MIP_wild_char = 'Z'

#aa_to_num = {'A':2,'M':1,'C':0}
#q=3
