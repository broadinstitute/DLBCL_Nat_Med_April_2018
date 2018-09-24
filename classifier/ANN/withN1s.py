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
        preds, reals, accuracy = trainNN.main(seed, useCOO, usen1s=True, originalReduced=originalReduced)
        if preds == prevPreds:
            print("seed issue")
        # print(preds)
        for i in range(0,len(preds)):
            allPreds.append(preds[i])
            allReals.append(reals[i])
        #allPreds.append(preds)
        #allReals.append(reals)
        accuracies.append(accuracy)
        prevPreds = preds

    return allPreds, allReals, accuracies

retvals = main(useCOO=False, numiter=3, originalReduced=True)

print("Accuracy: ", np.mean(retvals[2]))
pl.scatter(retvals[0], retvals[1])
pl.xlabel("Predicted Probability")
pl.ylabel("True Probability")
pl.title("Reduced Full Training Set N1s")
pl.savefig("reducedN1s.png")
pl.close()

print(len(retvals[0]))
print(pearsonr(retvals[0],retvals[1]))