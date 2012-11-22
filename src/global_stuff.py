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


base_folder = '../data/proteins/humvar/'
all_seqs_file = '../data/human-2011_12.seq'
neutral_mutations_file = '../data/humvar-2011_12.neutral.pph.input'
deleterious_mutations_file = '../data/humvar-2011_12.deleterious.pph.input'
#protein_list_file = '../data/hum_var_msa_completed'

protein_list_file = '/mnt/work/fultonw/deleterious/data/hum_var_real_msa_completed'

data_folder = '/mnt/work/fultonw/deleterious/data/'

completed_list_file = '../data/completed_list'

FILE_MANAGER_SIZE = 500
OBJ_MANAGER_SIZE = 500
#MUSCLE_PATH = '../muscle/muscle3.8.31_i86linux64'

MUSCLE_PATH = '/mnt/work/fultonw/muscle3.8.31_i86linux64'
DSSP_PATH = '/mnt/work/fultonw/active_site/dssp/dssp-2.0.4-linux-amd64'


BIN_FOLDER = '../data/bin/'
HOLDING_FOLDER = '../data/holding_folder/'
CACHE_MAX_SIZE = 500

to_reindex = True
recalculate = False

recalculate_nodewise_loss_f = True

metric_cutoffs = [1,2,3,4,5,6,7,8,9]

PROTEIN_BASE_FOLDER = '../data/proteins/humvar/'


RESULTS_FOLDER = '../new_results/testing_logreg/'
RESULTS_BASE_FOLDER = '../the_results/'

NACCESS_PATH = '/mnt/work/fultonw/active_site/Naccess/naccess'
NACCESS_FOLDER = '/mnt/work/fultonw/active_site/Naccess/'

#BLAST_PATH = '../blast/ncbi-blast-2.2.26+/bin/psiblast'
#BLASTDB_PATH = '../blast/ncbi-blast-2.2.26+/bin/nr/nr'

BLAST_PATH = '/mnt/work/fultonw/deleterious/blast/ncbi-blast-2.2.26+/bin/psiblast'
BLASTDB_PATH = '/mnt/work/fultonw/deleterious/blast/ncbi-blast-2.2.26+/bin/nr/nr'


CONSERVATION_FOLDER = '/home/fultonw/conservation_code/'

ORIG_CHAINS = '../catres_pdbs'
CSA_FILE = '../catres_sites'

success_file = 'success_catres.txt'
fail_file = 'fail_catres.txt'

proc_id = 0
