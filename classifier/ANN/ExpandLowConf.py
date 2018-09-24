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
        preds, reals, accuracy = trainNN.main(seed, useCOO, numExpand=5, lowExpanded=True, originalReduced=originalReduced)
        if preds == prevPreds:
            print("seed issue")
        # print(preds)
        allPreds.append(preds)
        allReals.append(reals)
        accuracies.append(accuracy)
        prevPreds = preds

    return allPreds, allReals, accuracies

retvals = main(useCOO=False, numiter=50, originalReduced=True)

print("Accuracy: ", np.mean(retvals[2]))
pl.scatter(retvals[0], retvals[1])
pl.xlabel("Predicted Probability")
pl.ylabel("True Probability")
pl.title("Reduced Full Training Set Low Expanded")
pl.savefig("reducedLowExpanded.png")
pl.close()