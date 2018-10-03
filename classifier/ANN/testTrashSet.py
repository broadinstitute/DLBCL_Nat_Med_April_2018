import torch
import torch.nn as nn
import NN
import glob
import pandas as pd
import numpy as np

netPaths = glob.glob('LatestModels/*')
nets = []

for currnet in netPaths:
    the_model = NN.Net(10,21,5)
    the_model.load_state_dict(torch.load(currnet))
    nets.append(the_model)

trashSet = pd.read_csv("DataTables/reducedDFTrash.tsv", sep = '\t')

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
filename = "NN_TrashSet_Confidences.tsv"
df.to_csv("WrittenFiles/"+filename, sep="\t")