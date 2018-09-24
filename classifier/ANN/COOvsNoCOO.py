import trainNN
import matplotlib.pyplot as pl
from scipy.stats.stats import pearsonr
import numpy as np

def main(useCOO, numiter, originalReduced):
    allPreds = []
    allReals = []
    accuracies = []
    prevPreds = None
    prevReals = None
    for i in range(0,numiter):
        seed = np.random.randint(low=0, high=999999)
        preds,reals,accuracy = trainNN.main(seed, useCOO, originalReduced=originalReduced)
        if preds == prevPreds:
            print("seed issue")
        #print(preds)
        allPreds.append(preds)
        allReals.append(reals)
        accuracies.append(accuracy)
        prevPreds = preds
        prevReals = reals

    return allPreds, allReals, accuracies

COOpoints = main(useCOO=True, numiter=50, originalReduced=False)
nonCOOpoints = main(useCOO=False, numiter=50, originalReduced=False)

print("Accuracy: ", np.mean(COOpoints[2]))
pl.scatter(COOpoints[0], COOpoints[1])
pl.xlabel("Predicted Probability")
pl.ylabel("True Probability")
pl.title("Reduced With COO")
pl.savefig("reducedwithCOO.png")
pl.close()

print("Accuracy: ", np.mean(nonCOOpoints[2]))
pl.scatter(nonCOOpoints[0], nonCOOpoints[1])
pl.xlabel("Predicted Probability")
pl.ylabel("True Probability")
pl.title("Reduced Minus COO")
pl.savefig("reducedminusCOO.png")
pl.close()