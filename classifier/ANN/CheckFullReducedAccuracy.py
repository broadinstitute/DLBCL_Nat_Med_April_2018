import trainNN
import matplotlib.pyplot as pl
from scipy.stats.stats import pearsonr
import numpy as np
import pandas as pd
import pylab

def main(useCOO, numiter, originalReduced):
    allPreds = []
    allReals = []
    accuracies = []
    sampleOrder = []
    prevPreds = None
    trainPreds, trainTruths = [],[]
    matches = []
    allPredLabels = []
    allActualLabels = []
    for i in range(0, numiter):
        seed = np.random.randint(low=0, high=999999)
        seed = 123
        preds, reals, accuracy, tmp1, tmp2, sampleNames, match, predLabels, actualLabels = trainNN.main(seed, useCOO, originalReduced=originalReduced)
        if preds == prevPreds and preds:
            print("seed issue")
        # print(preds)
        for i in range(0,len(preds)):
            allPreds.append(preds[i])
            allReals.append(reals[i])

        for i in range(0, len(sampleNames)):
            sampleOrder.append(sampleNames[i])
        accuracies.append(accuracy)
        prevPreds = preds
        for i in range(0,len(tmp1)):
            trainPreds.append(tmp1[i])
            trainTruths.append(tmp2[i])

        for i in range(0,len(predLabels)):
            allPredLabels.append(predLabels[i])
            allActualLabels.append(actualLabels[i])

        matches.append(match)

    return allPreds, allReals, accuracies, trainPreds, trainTruths, sampleOrder, matches, allPredLabels, allActualLabels

niter = 1
retvals = main(useCOO=False, numiter=niter, originalReduced=True)
colors = [item for sublist in retvals[6] for item in sublist]
allPredLabels = retvals[7]
allActualLabels = retvals[8]

print("Accuracy: ", np.mean(retvals[2]))
pl.scatter(retvals[1], retvals[0], c=colors)
pl.ylabel("Predicted Max Probability")
pl.xlabel("True Probability")
pl.title("Reduced Validation Set")
pylab.ylim([0,1])
pylab.xlim([0,1])
pl.savefig("reducedNoCOO.png")
pl.close()

#high confidence confmat
sPred = []
sAct = []
for i in range(0, len(allPredLabels)):
    if retvals[0][i] >= .90:
        sPred.append(allPredLabels[i])
        sAct.append(allActualLabels[i])
y_actu = pd.Series(sAct, name='Actual')
y_pred = pd.Series(sPred, name='Predicted')
confmat = pd.crosstab(y_pred,y_actu)
print(confmat)

pl.scatter(retvals[4], retvals[3])
pl.ylabel("Predicted Max Training Probability")
pl.xlabel("True Probability")
pl.title("Reduced Set Training Correlation")
pl.savefig("reducedNoCOOTraining.png")
pl.close()

#hash each person's absolute difference and then add to dataframe/sort
sampleList = retvals[5]
trainPreds = retvals[3]
trainTruths = retvals[4]
samplePredHash = {}
sampleTrueHash = {}
sampleDiffHash = {}
for i in range(0, len(sampleList)):
    currSample = sampleList[i]
    diff = np.abs(trainPreds[i]-trainTruths[i])
    if currSample not in sampleDiffHash:
        sampleDiffHash[currSample] = [diff]
        samplePredHash[currSample] = [trainPreds[i]]
        sampleTrueHash[currSample] = trainTruths[i]
    else:
        sampleDiffHash[currSample].append(diff)
        samplePredHash[currSample].append(trainPreds[i])

df = pd.DataFrame.from_dict(sampleDiffHash, orient='index')
df2 = pd.DataFrame.from_dict(samplePredHash, orient='index')
df3 = pd.DataFrame.from_dict(sampleTrueHash, orient='index')
dfAvg = df.copy()
dfAvg['AverageAbsoluteDiff'] = df.mean(numeric_only=True, axis=1)
dfAvg['stdOfDifferences'] = df.std(numeric_only=True, axis=1)
df2 = df2.reindex(dfAvg.index)
dfFull = pd.concat([df2,dfAvg], axis=1, sort=False)
colnames = []
for i in range(1,niter+1):
    colnames.append("max"+str(i))

for i in range(1,niter+1):
    colnames.append("diff"+str(i))

#colnames.append("TrueProbability")
colnames.append("AverageAbsoluteDiff")
colnames.append("stdOfDifferences")
colnames.append("TrueProbability")
df3 = df3.reindex(dfFull.index)
dfFull = pd.concat([dfFull,df3], axis=1, sort=False)

dfFull.columns = colnames

dfFull = dfFull.sort_values(by=['AverageAbsoluteDiff'], ascending=False)
dfFull = dfFull.iloc[:,::-1]


dfFull.to_csv("WrittenFiles/AbsoluteDifferences.csv")

print(len(retvals[0]), len(retvals[1]))
print(pearsonr(retvals[0],retvals[1]))