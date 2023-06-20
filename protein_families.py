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
#%% Lists of different protein family sequences. Each list has 200 seqs from Muscarinic acetylcholine receptor
#   Retinoid X receptor, and Retroviral VpR

MAR = ReadFastaFile('/home/research/Desktop/ARCC/cancerai/sequence(1).txt')
RXR = ReadFastaFile('/home/research/Desktop/ARCC/cancerai/sequence(2).txt')
RVP = ReadFastaFile('/home/research/Desktop/ARCC/cancerai/sequence(3).txt')

BCCT = ReadFastaFile('/home/research/Desktop/ARCC/cancerai/BCCT transporter family.txt')
FMR = ReadFastaFile('/home/research/Desktop/ARCC/cancerai/Flagellar M-ring protein(DNA).txt')
MET = ReadFastaFile('/home/research/Desktop/ARCC/cancerai/Metallothionein(DNA).txt')

Porin = ReadFastaFile('/home/research/Desktop/ARCC/cancerai/Porin_family.txt')

aquaPorin = ReadFastaFile("/home/research/Desktop/ARCC/cancerai/Aquaporin Z (all alphas).txt")
LamBPorin = ReadFastaFile("/home/research/Desktop/ARCC/cancerai/Porin, LamB-type (all betas).txt")
oprBPorin = ReadFastaFile("/home/research/Desktop/ARCC/cancerai/oprB porin (mix of betas and alphas).txt")


MAR_seq = MAR[:500]
RXR_seq = RXR[:500]
RVP_seq = RVP[:500]
BCCT_seq = BCCT[:500]
FMR_seq = FMR[:500]
MET_seq = MET[:500]
POR_seq = Porin[:500]

seq = MAR_seq + RXR_seq + RVP_seq
seq_2 = BCCT_seq + FMR_seq + MET_seq
seq_3 = seq + seq_2 + POR_seq + aquaPorin + LamBPorin + oprBPorin

#%% Deploy model on each protein family sequence:

import torch
import numpy as np
from torch.utils.data import DataLoader
from genslm import GenSLM, SequenceDataset

torch.manual_seed(1)

model = GenSLM("genslm_25M_patric", model_cache_dir="/home/research/Desktop/ARCC/cancerai")

device = torch.device("cuda:0")
model = model.to(device)
model.eval()

dataset = SequenceDataset(seq_3, model.seq_length, model.tokenizer)
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

'''
Below are two versions of PCA on the embeddings. 
'''

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
#from mpl_toolkits import mplot3d
#import mpld3

PCA = PCA(n_components=3)
pca_proj = PCA.fit(embeddings)
x = pca_proj.transform(embeddings)  

# build dataframe so that we can make the interactive plot:
df = pd.DataFrame({'x':x[:,0], 'y':x[:,1], 'z':x[:,2]})

fig = plt.figure(dpi=100)
ax = plt.axes(projection='3d')

ax.scatter3D(x[:500,0], x[:500,1], x[:500,2], label='MAR', s=6, color='blue')
ax.scatter3D(x[500:1000,0], x[500:1000,1], x[500:1000,2], label='RXR', s=4, color='red')
ax.scatter3D(x[1000:1500,0], x[1000:1500,1], x[1000:1500,2], label='RVP', s=2, color='green')
ax.scatter3D(x[1500:2000,0], x[1500:2000,1], x[1500:2000,2], label='BCCT', s=2, color='purple')
ax.scatter3D(x[2000:2500,0], x[2000:2500,1], x[2000:2500,2], label='FMR', s=2, color='orange')
ax.scatter3D(x[2500:3000,0], x[2500:3000,1], x[2500:3000,2], label='MET', s=2, color='black')


ax.set_title('Protein Families: PCA on all embeddings at once') 
ax.legend()

#%%
import plotly.graph_objects as go
from plotly.graph_objects import *
import plotly.express as px

fig = go.Figure(data=[go.Scatter3d(
    x = x[:500,0],
    y = x[:500,1],
    z = x[:500,2],
    mode ='markers',
    name = 'Muscarinic  Acetylcholine Receptor',
    marker = dict(
    size =1.5,
    color ='blue',  
        )),    
go.Scatter3d(
    x = x[500:1000,0],
    y = x[500:1000,1],
    z = x[500:1000,2],
    mode ='markers',
    name = 'Retroviral VpR',
    marker = dict(
    size = 1.5,
    color = 'red', 
        )),
go.Scatter3d(
    x = x[1000:1500,0],
    y = x[1000:1500,1],
    z = x[1000:1500,2],
    mode ='markers',
    name = 'Retinoid X Receptor',
    marker = dict(
    size = 1.5,
    color ='green',
        )),
go.Scatter3d(
    x = x[1500:2000,0],
    y = x[1500:2000,1],
    z = x[1500:2000,2],
    mode ='markers',
    name = 'BCCT Transporter',
    marker = dict(
    size = 1.5,
    color ='purple',    
        )),
go.Scatter3d(
    x = x[2000:2500,0],
    y = x[2000:2500,1],
    z = x[2000:2500,2],
    mode ='markers',
    name = 'Flagellar M-ring Protein (DNA)',
    marker = dict(
    size = 1.5,
    color ='orange',
        )),
go.Scatter3d(
    x = x[2500:3000,0],
    y = x[2500:3000,1],
    z = x[2500:3000,2],
    mode = 'markers',
    name = 'Metallothionein (DNA)',
    marker = dict(
    size = 1.5,
    color = 'black',
        )),
go.Scatter3d(
    x = x[3000:,0],
    y = x[3000:,1],
    z = x[3000:,2],
    mode = 'markers',
    name = 'Porin',
    marker = dict(
    size = 1.5,
    color = 'pink',
        )),
    ])

# tight layout
fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
fig.write_html("/home/research/Desktop/ARCC/cancerai/figure1.html") #Modifiy the html file
fig.show()








