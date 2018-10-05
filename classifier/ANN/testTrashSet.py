import torch
import torch.nn as nn
import NN
import glob
import pandas as pd
import numpy as np
import sys
import os

netPaths = glob.glob('LatestModels/*')
nets = []
fullFeatures = sys.argv[1] == "True"
if fullFeatures:
    nFeatures = 161
else:
    nFeatures = 21

for currnet in netPaths:
    the_model = NN.Net(10,nFeatures,5)
    the_model.load_state_dict(torch.load(currnet))
    nets.append(the_model)


if not fullFeatures:
    trashSet = pd.read_csv("DataTables/reducedDFTrash.tsv", sep = '\t')
else:
    trashSet = pd.read_csv("DataTables/trashV2.txt", sep = '\t')

trashSetVals = trashSet.values
trashSetTensors = []
for val in trashSetVals:
    valT = torch.tensor(np.asarray(val)).float().view(1, -1)
    trashSetTensors.append(valT)

outputs = []
for net in nets:
    for sample in trashSetTensors:
        outputs.append(net(sample).data[0].numpy())

df = pd.DataFrame(outputs)
filename = "NN_TrashSet_Confidences_features"+str(trashSet.shape[1])+".tsv"
df.to_csv("WrittenFiles/"+filename, sep="\t")