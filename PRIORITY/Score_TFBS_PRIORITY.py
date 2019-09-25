#!/usr/bin/env python
# coding: utf-8

# # Method 2 Site caller 
# (ALPHA version 1)
# 
# Author: Zachery Mielko
# Refactored by Tiffany Ho (Sept 16, 2019)
# 
# The script requires the following dependencies:
# - Python modules:
#     - pandas
#     - numpy
#     - biopython
# - Command line tools:
#     - Bedtools
# 
# 
# Method 2 takes in the following as **input**:
# 
# - Human Genome file (FASTA, must have an index file in the same directory, .fai
#     - You can get this from Samtools faidx
# - Alignment file from PRIORITY
# - Chip-seq peaks (BED file)
# 
# The script gives the following as **output**:
# - Bed file of centered sites (Centered_PRIORITY.bed)
# 
# The output as-is will just be the 1bp center, but you could extract the whole match, which is calculated. 
# 

# In[13]:


import pandas as pd
import numpy as np
import itertools
from Bio import SeqIO
import os, subprocess
import re
from Bio.Seq import reverse_complement
import pathlib

################## User defined parameters ###############
wkdir, TF = "/mnt/c/Users/th184/Box Sync/Gordan Lab/PRIORITY_score_caller", "Elk1" # run from Ubuntu hence /mnt/c
peak_coord_file = "Elk1DHSv2promoterHg19NoDac.bed" 
kmer_file = "Elk1_kmers_45_PRIORITY.txt"
genome_file = "/mnt/c/Users/th184/Box Sync/Probe Design/hg19.fa"

kmer_len = 7
overlap_req = 2
escore_cutoff = 0.45
# Prior assumptions about the kmer PRIORITY alignment
core = [7,12]  # 0-indexed
center_pos = 10


# In[2]:


io_folder = pathlib.Path(wkdir) / TF
os.chdir(io_folder)

# output - prob not here
called_TFBS_save = io_folder / f"called_TFBS_{str(escore_cutoff)}" 
centered_TFBS_save = io_folder / f"centered_TFBS_{str(escore_cutoff)}.bed"
centered_TFBS_verbose = io_folder / f"centered_TFBS_verbose_{str(escore_cutoff)}.bed"


# In[36]:


### kmer prep ###
# Some function requires kmers to be defined. SO this code defines the kmers from the PRIORITY alignment
kmers = pd.read_csv(io_folder/kmer_file , sep = '\t') # created an Unnamed column for no reason
kmers = kmers.loc[:, ~kmers.columns.str.contains('^Unnamed')] # get rid of the Unnamed column
iscore = []
position = []
kmer = []
for i in kmers['sequences']:
    if '.' in i[core[0]:core[1]+1]: # added +1
        iscore.append(False)
    else:
        iscore.append(True)
    kPos = re.search('(C|G|A|T)',i).start() # 0-indexed
    position.append(kPos) 
    kmer.append(i[kPos: kPos+kmer_len])
    #kmer.append(re.findall('[A-Z]{' + str(kmer_length) + '}',i)[0])
kmers['is_core'] = iscore
kmers['kPosition'] = position
kmers['kmer'] = kmer # already filtered to be >= threshold?!!
# kmers


# In[68]:


def fill_space(path):
    s = pathlib.PurePath(path).split()
    print(s)
    return "\ ".join((path).split())


# In[19]:


### getfasta prep ###
# Read CSV, filter for short seqs and duplicates, adjust coordinates
pk_coord = pd.read_csv(io_folder/peak_coord_file, sep = '\t', header = None, usecols=[0,1,2]) 
pk_coord = pk_coord[pk_coord[2]-pk_coord[1] > (kmer_len+overlap_req)].drop_duplicates() # filter out(?) short sequences the caller would have trouble with
pk_coord[1] = pk_coord[1] - 1 # Move back 1 to adjust for bedtools getfasta (0-indexed)

pk_coord_adj_file = peak_coord_file.split(".")[0] + "_0-indexed.bed"
pk_file = peak_coord_file.split(".")[0]+".fasta"
# need to figure the space issue out
pk_coord.to_csv(pk_coord_adj_file, sep = '\t', index = False, header = None) # output to getfasta
#subprocess.call(["bedtools", "getfasta", "-s", "-fi", fill_space(genome_file), "-bed", fill_space(io_folder)/pk_coord_adj_file, ">", fill_space(io_folder)/pk_file])
# subprocess.call(f'bedtools getfasta -s -fi {genome_file} -bed {io_folder/pk_coord_adj_file} > {io_folder/pk_file}')
# subprocess.call(f'bedtools getfasta -s -fi {fill_space(genome_file)} -bed {fill_space(pk_coord_adj_file)} > {fill_space(pk_file)}')
# subprocess.call(r"bedtools getfasta -s -fi /mnt/c/Users/th184/Box\ Sync/Probe\ Design/hg19.fa -bed /mnt/c/Users/th184/Box\ Sync/Gordan\ Lab/PRIORITY_score_caller/Elk1/Elk1DHSv2promoterHg19NoDac_0-indexed.bed > /mnt/c/Users/th184/Box\ Sync/Gordan\ Lab/PRIORITY_score_caller/Elk1/Elk1DHSv2promoterHg19NoDac.fasta")


# In[63]:



pk_file_test = "Elk1DHSv2promoterHg19NoDac_test.fasta"
def read_fasta(file):
    entry = []
    with open(pk_file_test, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            entry.append([record.id, str(record.seq).upper(), len(record.seq), str(reverse_complement(record.seq)).upper()])
    arr = np.array(entry)
    df = pd.DataFrame({'seq_id':arr[:, 0], 'seq':arr[:,1], 'seq_len':arr[:,2], 'rev_comp':arr[:,3]})
    return df
peak_df = read_fasta(pk_file_test)
peak_df


# In[67]:


def get_pos(seq, kmer_len, kmer_list):
    # return a list of pos at which kmer exists in the seq
    pos = []
    for i in range(len(seq) - kmer_len + 1):
        window = seq[i:i+kmer_len]
        if window in kmer_list.values:
            pos.append(i)
    return pos
            
def get_kmer_pos(peak_df, kmer_len, kmer_list):
    # add a colum to peak_df with pos of kmer within each seq
    peak_df["kmer_pos"] = peak_df.apply(lambda peak_df: get_pos(peak_df['seq'], kmer_len, kmer_list), axis = 1)
    return peak_df
    
get_kmer_pos(peak_df, kmer_len, kmers["kmer"])


# In[2]:


print("Reading Files...")
peaks = read_fasta(pk_file)
# Read fasta files
# PeaksWSeqs = pd.read_csv(pk_file, sep = '\t', header=None)
# PeaksWSeqs = DF_fasta(PeaksWSeqs)
# PeaksWSeqs['Length'] = PeaksWSeqs['Sequence'].apply(lambda x: len(x))
# Make a reverse complement of the sequences
# rc_PeaksWSeqs = PeaksWSeqs.copy(deep = True)
# rc_PeaksWSeqs['Sequence'] = rc_PeaksWSeqs['Sequence'].apply(lambda x: reverse_complement(x))
# turn the kmers series into a set for fast lookup
chosen_kmers = set(kmers['kmer'])
# Big function of the analysis
def PRIORITY_CALL(Sequences, orientation = '+'):
    Scored = Score_Blaster(x=Sequences, kmerSeqs = chosen_kmers) #Score sites
    Groups = Site_Caller(Scored) # Look for consecutive overlaps (consecutive in genome), groups are a unique value
    Called = Scored.copy(deep = True) # make a copy
    Called['Site_number'] = Groups # add overlap annotation to DF
    print("Finding Overlaps...")
    # Get overlaps of sites with more than 1 consecutive overlap
    values = Called.Site_number.value_counts() 
    Called = Called[Called.Site_number.isin(values.index[values.gt(overlap_req-1)])]
    Called['Position'] = Called['Position'].apply(lambda x: int(x))
    # Get the core kmers as a set. Core kmers are those that contain the core range
    core_kmers = set(kmers[kmers['is_core']]['kmer'])
    core_groups = []
    # for each group, only keep it if it has a core kmer
    for group in Called.groupby(by='Site_number'):
        if len(set(group[1].Seq).intersection(core_kmers)) > 0:
            core_groups.append(group[0])
    # Get the core groups as a set, filter by core groups
    core_groups = set(core_groups)
    Called_core = Called[Called['Site_number'].isin(core_groups)]
    # Add kPosition to the DF
    Called2 = pd.merge(Called_core, kmers[['kmer','kPosition']], left_on = 'Seq', right_on='kmer')
    Called2 = Called2.sort_values(by=['Site_number', 'Position'])
    # Check that the kPosition is consecutive
    def checkConsecutive(l): 
        return sorted(l) == l 
    con_groups = []
    for group in Called2.groupby(by='Site_number'):
        if checkConsecutive(list(group[1].kPosition)):
            con_groups.append(group[0])
    con_groups = set(con_groups)
    Called2 = Called2[Called2['Site_number'].isin(con_groups)]
    Called2 = pd.merge(Called2, PeaksWSeqs[['Position','Length']], left_on = 'Peak_Seq', right_on = 'Position')
    Called3 = Chrom_Splitter(Called2)
    final = Genomic_Adjuster(Called3, orient=orientation)
    final = final.drop_duplicates()
    return(final)
# Run the function for positive and negative strands
Positive = PRIORITY_CALL(PeaksWSeqs)
Negative = PRIORITY_CALL(rc_PeaksWSeqs, orientation= '-')
# Merge the results
total = pd.concat([Positive,Negative])
total['Start'] = total['Center']
total['End'] = total['Center']
total = total[['Chromosome', 'Start', 'End', 'Orient']]
total.to_csv(f"{io_folder}/Centered_PRIORITY.bed", sep = '\t', header = None, index = False)


# In[75]:


##### ORIGINAL CODE #####

def DF_fasta(x):
    pos = x[x[0].str.startswith('>')]
    pos = pos.reset_index(drop=True)
    seqs = x[~x[0].str.startswith('>')]
    seqs = seqs.reset_index(drop=True)
    new = pd.DataFrame({'Position':pos[0],'Sequence':seqs[0]})
    return(new)
def read_fasta(file): # never called
  records = []
  sequence = []
  with open(file, "rU") as handle:
      for record in SeqIO.parse(handle, "fasta"):
        records.append(record.id)
        sequence.append(str(record.seq))
  Whole_seqs = pd.DataFrame({'Name':records,'Sequence':sequence})
  return Whole_seqs

def kmer_match(seq, kmer_len, kmer_list, Item):
    """kmer_match takes a string and scans along it to match the score with E-scores"""
    # return all the positions within the seq at which the kmer is within the kmer_list (which >= threshold)
    seq_list = []
    pos_list = []
    for i in range(len(String) - (window-1)):
        seq = String[i:i+window]
        pos = i
        if seq in kmer_list: 
## DOES NOT WORK! using in check with DataFrame and Series checks whether the val is contained in the Index!!!
            seq_list.append(seq)
            pos_list.append(pos)
    result = pd.DataFrame({'Seq':seq_list,'Position':pos_list,'Peak_Seq':Item})
    return(result)

def Score_Blaster(x, kmerSeqs,Thresh=Threshold, kmer_length = kmer_length):
    """Score_Blaster applies kmer_match to a dataframe. Only gives scores above a threshold, zip version"""
    # apply kmer_match to ALL seq
## threshold NEVER USED
    DataFrames = []
    print('Total length: ' + str(len(x)))
    for index,row in enumerate(zip(x['Sequence'], x['Position'])):
      # Counter to keep track of progress
      if int(index)%5000 == 0:
        print(index)
      # kmer match function 
      DataFrames.append(kmer_match(row[0], kmer_length, kmerSeqs,  row[1])) #SUPER SLOOOOOW
    result = pd.concat(DataFrames)
    return(result)

def get_kmer_pos_OLD(peak_df, kmer_len, kmer_list): # add Item (aka pos in peak_df) as needed
    # meant to combine kmer_match and score_blaster
    # decide not to use b/c appraoch  unnecessarily convoluted
    # replaced with get_kmer_pos()
    '''
    for each sequence within the peak_df, get the pos of kmers with E-score >= cutoff
    '''
    seq_list, pos_list = [], []
    for row in zip(peak_df['seq'], peak_df['seq_id']):
        seq, pos = row[0], row[1]
        #print(seq, pos)
        for i in range(len(seq) - kmer_len + 1):
            window = seq[i:i+kmer_len]
            
            if window in kmer_list.values:
                seq_list.append(window)
                pos_list.append(i) # maybe add "pos" as well?
    df = pd.DataFrame(list(zip(seq_list, pos_list)), columns =['Peak_Seq', 'Position']) # ['seq', 'pos']
    return df

def Site_Caller(x):
  """ Siter caller looks for consecutive overlaps """
  Groups = []
  PrevPos = x['Position'][0] -1
  PrevSeq = x['Peak_Seq'][0]
  SeqNumber = 0
  for idx,row in enumerate(zip(x['Peak_Seq'],x['Position'])):
    # When a new Sequence is being read
    if str(row[0]) != str(PrevSeq):
      PrevSeq = str(row[0])
      PrevPos = row[1]
      SeqNumber = SeqNumber + 1
      Groups.append(SeqNumber)
    # When you have a consecutive position
    elif row[1] - PrevPos == 1:
      Groups.append(SeqNumber)
      PrevPos = row[1]
    # When you have a non-consecutive position
    elif row[1] - PrevPos != 1:
      PrevPos = row[1]
      SeqNumber = SeqNumber + 1
      Groups.append(SeqNumber)
  return(Groups)


def Chrom_Splitter(x):
    """Chrom_Splitter takes the concatinated names given in fasta outputs from bedtool's getfasta and turns them into bed compatible columns"""
    Chrom = []
    Start = []
    End = []
    for i in x.Peak_Seq:
        i = i[1:-2]
        cr = i.split(':')
        pos = cr[1].split('-')
        Chrom.append(cr[0])
        Start.append(int(pos[0]))
        End.append(int(pos[1]))
    x['Chromosome'] = Chrom
    x['Start'] = Start
    x['End'] = End
    return(x)

def Genomic_Adjuster(x,orient, kmer_length=kmer_length):
    g = x.groupby(by=['Peak_Seq', 'Site_number'])
    Site_Start = []
    Site_End = []
    Chrom = []
    Center = []
    for name, group in g:
        if orient == '+':
            group_start = int(min(group.Position_x)) + int(group.Start.unique())
            group_end = int(max(group.Position_x)) + int(group.Start.unique()) + kmer_length
            group_center = (center_pos - int(min(group.kPosition))) + int(group.Start.unique()) + int(min(group.Position_x))
            Site_Start.append(group_start)
            Site_End.append(group_end)
            Center.append(group_center)
            Chrom.append(group.Chromosome.unique()[0])
        elif orient == '-':
            slength = int(group.Length.unique())
            group_end = int(slength - min(group.Position_x)) + int(group.Start.unique())
            group_start = int(slength - max(group.Position_x) - kmer_length) + int(group.Start.unique()) 
            group_center = (slength - (center_pos - int(min(group.kPosition)) + int(min(group.Position_x)))) + int(group.Start.unique()) +1
            Site_Start.append(group_start)
            Site_End.append(group_end)
            Center.append(group_center)
            Chrom.append(group.Chromosome.unique()[0])
    result = pd.DataFrame({'Chromosome':Chrom, "Start":Site_Start,
             'End':Site_End, 'Center':Center, 'Orient':orient})
    print(f"Orientation selected is: {orient}")
    return(result)

