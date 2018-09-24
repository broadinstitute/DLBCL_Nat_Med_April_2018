import trainNN_kfold
import matplotlib.pyplot as pl
from scipy.stats.stats import pearsonr
import numpy as np
import pandas as pd
import pylab
import torch
import sys

useSV = sys.argv[1] == "True"
useCNA = sys.argv[2] == "True"
downSample = False

def main(seed, folds):
    np.random.seed(seed)
    if((len(sys.argv) != 3)):
        print("provide two True/False arguments: first for useSV and second for useCNA")
        exit()
    nets, validationFrames, validationLabels = trainNN_kfold.main(seed, folds,
                                                                  original=True,
                                                                  downSample=downSample,
                                                                  useSV=useSV,
                                                                  useCNA=useCNA)

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
            idx = np.argmax(validationLabels[k][i])
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

    return predictions, trueLabels, probabilities, match, sampleNames

folds = 4
#seed = 1234
seed = np.random.randint(0,99999999)
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
if(downSample):
    filename = str(folds)+"foldNN_DS.csv"
else:
    if not useSV and not useCNA:
        filename = str(folds) + "foldNN_noSV_noCNA.csv"
    elif not useCNA:
        filename = str(folds) + "foldNN_noCNA.csv"
    elif not useSV:
        filename = str(folds) + "foldNN_noSV.csv"
    else:
        filename = str(folds) + "foldNN.csv"
df.to_csv("WrittenFiles/"+filename)
df.to_csv("/Users/twood/RProjects/LymphomaClean/DataTables/"+filename)

