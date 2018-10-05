import trainNN_kfold
import matplotlib.pyplot as pl
from scipy.stats.stats import pearsonr
import numpy as np
import pandas as pd
import pylab
import torch
import sys
import os
import glob

print(sys.argv)
useSV = sys.argv[1] == "True"
useCNA = sys.argv[2] == "True"
GD = sys.argv[3] == "True"
useProbs = sys.argv[4] == "True"
fixedValidationSets = sys.argv[5] == "True"
fullFeatures = sys.argv[6] == "True"
downSample = False

def main(seed, folds):
    np.random.seed(seed)
    if((len(sys.argv) != 7)):
        print("provide five True/False arguments: first for useSV, second for useCNA, third for genome doubling, \n"
              "fourth for probabilistic labels, fifth for fixed validation sets")
        exit()
    nets, validationFrames, validationLabels = trainNN_kfold.main(seed, folds,
                                                                  original=True,
                                                                  downSample=downSample,
                                                                  useSV=useSV,
                                                                  useCNA=useCNA,
                                                                  GD = GD,
                                                                  useProbs=useProbs,
                                                                  fixedValidationSets=fixedValidationSets,
                                                                  fullFeatures=fullFeatures)
    print(validationFrames[1].shape)
    validationValues = []
    for k in range(0, folds):
        vals = validationFrames[k].values
        tensorVals = []
        for val in vals:
            valT = torch.tensor(np.asarray(val)).float().view(1, -1)
            tensorVals.append(valT)
        validationValues.append(tensorVals)

    predictions = []
    probabilities = []
    for k in range(0, folds):
        net = nets[k]
        preds = []
        probs = []
        for i in range(0, len(validationValues[k])):
            output = net(validationValues[k][i]).data
            idx = torch.max(output, 1)
            prob = idx[0].item()
            idx = idx[1].item()
            preds.append(idx)
            probs.append(prob)
        predictions.append(preds)
        probabilities.append(probs)

    trueLabels = [[] for k in range(folds)]
    match = [[] for k in range(folds)]
    for k in range(0, folds):
        for i in range(0, len(validationLabels[k])):
            if useProbs:
                idx = np.argmax(validationLabels[k][i])
            else:
                idx = validationLabels[k][i][0]
            trueLabels[k].append(idx)
            if(idx == predictions[k][i]):
                match[k].append(True)
            else:
                match[k].append(False)

    sampleNames = []
    for k in range(folds):
        currFrame = validationFrames[k]
        for idx in currFrame.index:
            sampleNames.append(idx)

    files = glob.glob('LatestModels/*')
    for f in files:
        os.remove(f)
    i = 0
    nFeatures = validationFrames[1].shape[1]
    for currnet in nets:
        netName = "Net_seed"+str(seed)+"_Fold"+str(i)+"_nFeatures"+str(nFeatures)+".NN"
        torch.save(currnet.state_dict(), "Models/" + netName)
        torch.save(currnet.state_dict(), "LatestModels/" + netName)
        i += 1

    return predictions, trueLabels, probabilities, match, sampleNames

folds = 4
seed = 44825199
#seed = np.random.randint(0,99999999)
predictions, trueLabels, probabilities, matches, sampleNames = main(seed, folds)
predictionsFlat = [item for sublist in predictions for item in sublist]
trueLabelsFlat = [item for sublist in trueLabels for item in sublist]
probabilitiesFlat = [item for sublist in probabilities for item in sublist]
matchesFlat = [item for sublist in matches for item in sublist]


#add 1 to fix cluster labels
predictionsFlat = [i+1 for i in predictionsFlat]
trueLabelsFlat = [i+1 for i in trueLabelsFlat]
df = pd.DataFrame({
    'Sample' : sampleNames,
    'maxVal' : probabilitiesFlat,
    'correct' : matchesFlat,
    'predictedCluster' : predictionsFlat,
    'trueCluster' : trueLabelsFlat
})
df = df.sort_values(by='maxVal')
df = df.reset_index(drop=True)
filename = str(folds)+"foldNN"
if(downSample):
    filename = filename+"_DS"
if not useSV:
    filename = filename+"_noSV"
if not useCNA:
    filename = filename+"_noCNA"
if GD:
    filename = filename+"_GD"
if not useProbs:
    filename = filename+"_BinaryLabels"
if fixedValidationSets:
    filename = filename+"_FixedSets"
if fullFeatures:
    filename = filename+"_161F"

filename = filename+".csv"
df.to_csv("WrittenFiles/"+filename)
