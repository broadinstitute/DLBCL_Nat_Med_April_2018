import trainNN
import matplotlib.pyplot as pl
from scipy.stats.stats import pearsonr
import numpy as np

def main(useCOO, numiter, originalReduced):
    allPreds = []
    allReals = []
    accuracies = []
    prevPreds = None
    for i in range(0, numiter):
        seed = np.random.randint(low=0, high=999999)
        preds, reals, accuracy, finalPredictions, finalTruths, sampleNames = trainNN.main(seed, useCOO, numExpand=3, highExpanded=True, originalReduced=originalReduced)
        if preds == prevPreds:
            print("seed issue")
        for i in range(0,len(preds)):
            allPreds.append(preds[i])
            allReals.append(reals[i])
        accuracies.append(accuracy)
        prevPreds = preds
        print(accuracies)

    return allPreds, allReals, accuracies

retvals = main(useCOO=False, numiter=1, originalReduced=True)

print("Accuracy: ", np.mean(retvals[2]))
pl.scatter(retvals[1], retvals[0])
pl.ylabel("Predicted Probability")
pl.xlabel("True Probability")
pl.title("Reduced Full Training Set High Expanded")
pl.savefig("reducedHighExpanded.png")
pl.close()

print(pearsonr(retvals[0],retvals[1]))