import trainNNFull
import matplotlib.pyplot as pl
from scipy.stats.stats import pearsonr
import numpy as np

def main(useCOO, numiter, originalReduced):
    allPreds = []
    allReals = []
    accuracies = []
    prevPreds = None
    trainPreds, trainTruths = None,None
    for i in range(0, numiter):
        seed = np.random.randint(low=0, high=999999)
        #fixed seed
        seed = 144533
        preds, reals, accuracy, tmp1, tmp2 = trainNNFull.main(seed)
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
        trainPreds = tmp1
        trainTruths = tmp2

    return allPreds, allReals, accuracies, trainPreds, trainTruths

retvals = main(useCOO=False, numiter=1, originalReduced=True)

print("Accuracy: ", np.mean(retvals[2]))
pl.scatter(retvals[1], retvals[0])
pl.ylabel("Predicted Probability")
pl.xlabel("True Probability")
pl.title("FullDF Full Training Set")
pl.savefig("fullN1s.png")
pl.close()

pl.scatter(retvals[4], retvals[3])
pl.ylabel("Predicted Training Probability")
pl.xlabel("True Probability")
pl.title("FullDF Full Training Correlation")
pl.savefig("fullN1sTraining.png")
pl.close()

print(len(retvals[0]), len(retvals[1]))
print(pearsonr(retvals[0],retvals[1]))