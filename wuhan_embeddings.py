# This is a function that will convert a fasta file into a list

def ReadFastaFile(filename):
    fileObj = open(filename, 'r')
    sequences = []
    seqFragments = []
    for line in fileObj:
        if line.startswith('>'):
            if seqFragments:
                sequence = ''.join(seqFragments)
                sequences.append(sequence)
            seqFragments = []
        else:
            seq = line.rstrip()
            seqFragments.append(seq)
    if seqFragments:
        sequence = ''.join(seqFragments)
        sequences.append(sequence)
    fileObj.close()
    return sequences

#%%

# we are artificially creating mutations in two strains: Wuhan and Guangdong. We use Guangdong as a proxy to really see
# the behavior of the Wuhan strain in the embedding space.

Wuhan_seq = ReadFastaFile('/home/research/Desktop/ARCC/BVBRC_genome_sequence.fasta')

#%%

# creating a dictionary for amino acids:
import pandas as pd

df = pd.read_csv('/home/research/Desktop/ARCC/codons.tsv', header=None, delimiter=' ')    

AA_dict = dict(zip(df[0], df[1]))

#%%

from itertools import product

Wstring = Wuhan_seq[0].upper()

# changed the DNA sequence into codons to identify amino acids:
Wseq = Wstring.replace("T", "U")
    
bp = ['G', 'A', 'U', 'C']

codons = [''.join(comb)[:3] for comb in product(bp, repeat=len(bp))]
codons = list(set(codons))

# want to remove the stop codons so we can just focus on syn/nonsyn mutations
codons.remove('AGA')
codons.remove('AGG')
codons.remove('UGA')
#print(codons)

# want to break apart sequence so that we can substitute codons
Wseq_ = [Wseq[i:i+3] for i in range(0, len(Wseq), 3)]
# try w/ random probability, replace codons. Select 4 codons that include syn/nonsyn mutations.
# then generate ~100 strains, try pca, and interpret results...

# convert bp seq to AA seq:
    
AA_Wseq = [AA_dict.get(x,x) for x in Wseq_]

AA_Wseq_ = [AA_dict.get(x,x) for x in Wseq_]

# generate random mutations
import random

def rand_mut_index_generator(num_substitutions, num_mut_strs):
    '''
    Parameters
    ----------
    num_substitutions : int
        Number of mutations you would like in the DNA sequence.
    num_mut_strs : int
        Number of mutated DNA strains you want.

    Returns
    -------
    indexes : A list where each element is a list that contains index positions
              based on how many substitution mutations you want.
    '''
    indexes = []
    for i in range(0,num_mut_strs):
        rand_list = [random.randint(0, len(Wseq_)) for i in range(0,num_substitutions+1)]
        rand_list.sort()
        indexes.append(rand_list)
    return indexes

num_substitutions = 200
num_mut_strs = 300

indexes = rand_mut_index_generator(num_substitutions, num_mut_strs)
#%%  
    
def mutated_seqs_generator(og_seq, num_mut_strs, num_substitutions, indexes, codon):
    '''
    Parameters
    ----------
    og_seq : List
        A list of the original DNA sequence but each element is a codon string.
    num_mut_strs : int
        Number of mutated strains you want.
    num_substitutions : int
        Number of mutations you would like in the DNA sequence.
    indexes : List of Lists
        A list where each element is a list that contains the number of 
        index positions that you want to perform substitution.
    codon : str
        A codon in the form of a str. Ex: 'UUU'

    Returns
    -------
    mut_seqs : A list of mutated strains derived from the orignal sequence.
    '''
    rand_mut_strs = []
    for i in range(0, num_mut_strs):
        og_seq_ = [og_seq[k:k+3] for k in range(0, len(og_seq), 3)]
        for j in range(0,num_substitutions):
            og_seq_[indexes[i][j]] = codon
        rand_mut_strs.append(og_seq_)
        
    mut_seqs = []
    for lists in rand_mut_strs:
        seq = ''.join(lists)
        mut_seqs.append(seq)  
    return mut_seqs

mutated_strains = mutated_seqs_generator(Wseq, num_mut_strs, num_substitutions, indexes, codon='GUC')

def codon_to_AA(codon_seq, AA_dict):
    '''
    Parameters
    ----------
    seq : List
        The DNA sequence split up into elements of codons in the form of strs
    AA_dict : Dictionary
        A dictionary that contains the codon corresponing to an amino acid

    Returns
    -------
    AA_seq : A list that contains the converted amino acids
    '''
    seq_ = [codon_seq[i:i+3] for i in range(0, len(codon_seq), 3)]
    AA_seq = [AA_dict.get(x,x) for x in seq_]
    return AA_seq
#%%

mut_str_seq_ = [mutated_strains[0][i:i+3] for i in range(0, len(mutated_strains[0]), 3)]
mut_str_AA_seq_1 = codon_to_AA(mutated_strains[0], AA_dict)
synnonsyn_index = [x for x in range(0,len(mut_str_seq_))]


codon_AA_tuple = list(zip(synnonsyn_index, mut_str_seq_, mut_str_AA_seq_1))

W_codon_AA_tuple = list(zip(synnonsyn_index, Wseq_, AA_Wseq_))
#%%
synnonsyn_list = []

for i in range(0, len(codon_AA_tuple)):
    if codon_AA_tuple[i][2] != W_codon_AA_tuple[i][2]:
        synnonsyn_list.append('nonsynonymous')
    else:
        synnonsyn_list.append('synonymous')

synnonsyn_tuple = list(zip(synnonsyn_index, synnonsyn_list))


syn_mutations = []
nonsyn_mutations_index = []

# note, that from this block below, it appears that some of the information on the remaining syn mutations is lost
# which i think is due to the fact that maybe some of the substitution mutations overwrote themselves...

for i in range(0, len(codon_AA_tuple)):
    if codon_AA_tuple[i][2] != W_codon_AA_tuple[i][2]:
        nonsyn_mutations_index.append(codon_AA_tuple[i][0])
    elif codon_AA_tuple[i][2] == W_codon_AA_tuple[i][2] and codon_AA_tuple[i][1] != W_codon_AA_tuple[i][1]:
        syn_mutations.append(codon_AA_tuple[i])

new_Wseq_ = [Wseq[i:i+3] for i in range(0, len(Wseq), 3)]

for i in nonsyn_mutations_index:
    new_Wseq_[i] = 'GUC'

new_Wseq_=''.join(new_Wseq_)




#%% plotting different codon substitution strain embeddings:
    
mutated_strains = mutated_seqs_generator(Wseq, num_mut_strs, num_substitutions, indexes, codon='GUC')
mutated_strains1 = mutated_seqs_generator(Wseq, num_mut_strs, num_substitutions, indexes, codon='UUU')
mutated_strains2 = mutated_seqs_generator(Wseq, num_mut_strs, num_substitutions, indexes, codon='UAA')

embedding_list = mutated_strains + mutated_strains1 + mutated_strains2

#%%

import torch
import numpy as np
from torch.utils.data import DataLoader
from genslm import GenSLM, SequenceDataset

torch.manual_seed(1)

model = GenSLM("genslm_25M_patric", model_cache_dir="/home/research/Desktop/ARCC/cancerai")

device = torch.device("cuda:0")
model = model.to(device)
model.eval()

dataset = SequenceDataset(embedding_list, model.seq_length, model.tokenizer)
dataloader = DataLoader(dataset, batch_size = 6, shuffle=False)

# Compute averaged-embeddings for each input sequence
embeddings = []

with torch.no_grad():
    for batch in dataloader:
        outputs = model(batch["input_ids"].to(device), batch["attention_mask"].to(device), output_hidden_states=True)
        # outputs.hidden_states shape: (layers, batch_size, sequence_length, hidden_size)
        emb = outputs.hidden_states[0].detach().cpu().numpy()
        # Compute average over sequence length
        emb = np.mean(emb, axis=1)
        embeddings.append(emb)
    torch.cuda.empty_cache()
    
# Concatenate embeddings into an array of shape (num_sequences, hidden_size)
embeddings = np.concatenate(embeddings)
embeddings.shape

#%%

# PCA on the embeddings

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

PCA = PCA(n_components=2)
pca_proj = PCA.fit(embeddings)
x = pca_proj.transform(embeddings)  

fig, ax = plt.subplots(figsize = (8,8))

ax.scatter(x[:300,0], x[:300,1], label='GUC', s=6, color='blue')
ax.scatter(x[300:600,0], x[300:600,1], label='UUU', s=4, color='red')
ax.scatter(x[600:,0], x[600:,1], label='UAA', s=2, color='green')

ax.set_title('Wuhan') 
ax.legend()















