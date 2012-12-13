'''
Created on Mar 9, 2012

@author: glareprotector
'''
#from manager import *
import string
#import Bio.PDB
import csv
#import constants
import string
import re
import math

import pdb

def get_param():
    import param
    p = param.param({'ev':1e-10, 'protein_list_file':'mf_done', 'uniprot_id':'Q8WXA2', 'avg_deg':3, 'n_cutoff':0, 'f_cutoff':15, 'which_msa':1, 'which_weight':1, 'which_dist':1, 'pseudo_c':1})
    return p

#home = '/home/fw27/d/deleterious/'
home = '/mnt/work/fultonw/deleterious/'

polyphen_msa_directory = home+'data/polyphen-2.2.2/precomputed/alignments/'


lock_folder = home + 'lock_folder/'

#base_folder = '../data/proteins/humvar_from_orchestra/'
base_folder = '/mnt/work/fultonw/deleterious/data/proteins/humvar/'
#base_folder = home

all_seqs_file = '../data/human-2011_12.seq'
neutral_mutations_file = '../data/humvar-2011_12.neutral.pph.input'
deleterious_mutations_file = '../data/humvar-2011_12.deleterious.pph.input'
#protein_list_file = '../data/hum_var_msa_completed'


cosmic_gene_list =home+ 'data/cosmic_genes'
cosmic_raw_data_folder = home+'data/fasta/'

protein_list_file = home+'data/humvar_list'

data_folder = home+'data/'

completed_list_file = '../data/completed_list'

FILE_MANAGER_SIZE = 500
OBJ_MANAGER_SIZE = 500
#MUSCLE_PATH = '../muscle/muscle3.8.31_i86linux64'

MUSCLE_PATH = home+'muscle/muscle3.8.31_i86linux64'
DSSP_PATH = '/mnt/work/fultonw/active_site/dssp/dssp-2.0.4-linux-amd64'


BIN_FOLDER = '../data/bin/'
HOLDING_FOLDER = '../data/holding_folder/'
CACHE_MAX_SIZE = 50

to_reindex = True
recalculate = False

recalculate_nodewise_loss_f = True

metric_cutoffs = [1,2,3,4,5,6,7,8,9]

#aa_to_num = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'J':9,'K':10,'L':11,'M':12,'N':13,'O':14,'P':15,'Q':16,'R':17,'S':18,'T':19,'U':20,'V':21,'W':22,'X':23,'Y':24,'Z':25,'-':26}

aa_to_num = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'-':20}

q = 21



#PROTEIN_BASE_FOLDER = home+'data/proteins/humvar_from_orchestra/'
PROTEIN_BASE_FOLDER = base_folder

COSMIC_BASE_FOLDER = home+'data/proteins/cosmic/'

dropbox_folder = '~/Dropbox/deleterious/'

RESULTS_FOLDER = '../new_results/testing_logreg/'
RESULTS_BASE_FOLDER = '../the_results/'

NACCESS_PATH = '/mnt/work/fultonw/active_site/Naccess/naccess'
NACCESS_FOLDER = '/mnt/work/fultonw/active_site/Naccess/'

#BLAST_PATH = '../blast/ncbi-blast-2.2.26+/bin/psiblast'
#BLASTDB_PATH = '../blast/ncbi-blast-2.2.26+/bin/nr/nr'



BLAST_PATH = '/mnt/work/fultonw/deleterious/blast/ncbi-blast-2.2.26+/bin/psiblast'
#BLASTDB_PATH = '/mnt/work/fultonw/nr/'
BLASTDB_PATH = 'nr/nr'

CONSERVATION_FOLDER = '/home/fultonw/conservation_code/'

ORIG_CHAINS = '../catres_pdbs'
CSA_FILE = '../catres_sites'

success_file = 'success_catres.txt'
fail_file = 'fail_catres.txt'

proc_id = 0
