import NN
import GenerateTestTrainSets as gt
import numpy as np
import torch
import torch.nn as nn
import pandas as pd
import matplotlib.pyplot as pl
from scipy.stats.stats import pearsonr

def main(seed, useCOO, originalReduced, usen1s=False, numExpand=0, highExpanded=False, lowExpanded=False, verbose=False):
    useExpandedSets = False
    if useExpandedSets:
        trainSet, validationSet, trainLabels, validationLabels = gt.generateExpandedSets(expandSize=1)
        validationLabels = np.array(validationLabels) - 1
    elif lowExpanded:
        trainSet, validationSet, trainLabels, validationLabels = gt.generateLowExpandedDF(seed=seed, numExpand=numExpand,
                                                                                originalReduced=originalReduced,
                                                                                useProbs=True, COO=False)
    elif highExpanded:
        trainSet, validationSet, trainLabels, validationLabels = gt.generateHighExpandedDF(seed=seed,
                                                                                          numExpand=numExpand,
                                                                                          originalReduced=originalReduced,
                                                                                          useProbs=True, COO=False)
    else:
        trainSet, validationSet, trainLabels, validationLabels = gt.generateDFs(seed = seed,originalReduced=originalReduced,usen1s=usen1s,useProbs=True, COO=useCOO)



    torch.manual_seed(333)

    net = NN.Net(10,trainSet.shape[1],5)
    print("Training size: " + str(trainSet.shape[0]))
    trainValues = trainSet.values
    validationValues = validationSet.values
    sampleNames = trainSet.index.values


    # lr = 0.01
    # num_epoch = 200

    lr = 0.1
    num_epoch = 100

    MSE = True

    if MSE:
        criterion = nn.MSELoss()
    else:
        criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.SGD(net.parameters(), lr=lr, weight_decay=.001)
    #optimizer = torch.optim.ASGD(net.parameters())

    trainInputsArr = []
    trainLabelsArr = []

    for i in range(0,len(trainValues)):
        if not MSE:
            toAppend = torch.tensor([trainValues[i]], requires_grad=True).float()
            toAppend2 = torch.tensor([int(trainLabels[i]-1)]).long()
        else:
            toAppend = torch.tensor([trainValues[i]], requires_grad=True).float()
            toAppend2 = torch.tensor([trainLabels[i]]).float()
        trainLabelsArr.append(toAppend2)
        trainInputsArr.append(toAppend)

    finalPredictions = []
    finalTruths = []
    for epoch in range(0, num_epoch):
        avgloss = 0
        currPreds = []
        currTruths = []
        for i in range(0, len(trainInputsArr)):
            X = trainInputsArr[i]
            Y = trainLabelsArr[i]
            #feedforward - backprop
            optimizer.zero_grad()
            out = net.forward(X)
            loss = criterion(out, Y)
            avgloss += loss.item()
            loss.backward()
            optimizer.step()
            currp = torch.max(Y, 1)[0].item()
            #idx = torch.max(Y, 1)[1].item()
            val = np.max(out.data[0].numpy())
            #val = val[idx]
            currPreds.append(val)
            currTruths.append(currp)

        avgloss = avgloss/len(trainInputsArr)
        finalPredictions = currPreds
        finalTruths = currTruths

        toprint = True
        if(toprint):
            if (epoch) % (num_epoch/10) == 0:
                print ('Epoch [%d/%d] Loss: %.4f' % (epoch, num_epoch, avgloss))


    validationIns = []
    validationLabs = []
    for i in range(0,len(validationValues)):
        toAppend = torch.tensor(np.asarray((validationValues[i]))).float().view(1,-1)
        validationIns.append(toAppend)

    preds = []
    probs = []
    for i in range(0,len(validationIns)):
        output = net(validationIns[i]).data
        idx = torch.max(output,1)
        prob = idx[0].item()
        idx = idx[1].item()
        preds.append(idx)
        probs.append(prob)

    if verbose:
        if not MSE:
            print(preds)
            print(validationLabels)

    correct = 0
    match = []
    for i in range(0, len(preds)):
        if(preds[i] == np.argmax(validationLabels[i])):
            correct = correct+1
            match.append("b")
        else:
            match.append("r")

    if verbose:
        for i in range(0,len(probs)):
            if(preds[i] != np.argmax(validationLabels[i])):
                print(preds[i],validationLabels[i], probs[i])

    if preds:
        accuracy = correct/len(preds)
    else:
        accuracy = 0

    if verbose:
        print("Accuracy: "+str(accuracy))

    sPred = pd.Series(preds, name="Prediction")
    correctLabels = [np.argmax(x) for x in validationLabels]
    sAct = pd.Series(correctLabels, name="Actual")
    confmat = pd.crosstab(sPred,sAct)
    if verbose:
        print(confmat)

    correctProbs = [max(x) for x in validationLabels]
    if verbose:
        print(probs)
        print(correctProbs)

    diffList = []
    for i in range(0,len(probs)):
        diffList.append(np.abs(probs[i] - correctProbs[i]))

    if verbose:
        print(diffList)
        print(np.mean(probs))
        print(np.mean(correctProbs))
        print(np.mean(diffList))

        print(pearsonr(probs,correctProbs))

    # pl.scatter(probs,correctProbs)
    # pl.xlabel("Predicted Probability")
    # pl.ylabel("True Probability")
    # pl.show()
    return probs,correctProbs,accuracy,finalPredictions,finalTruths,sampleNames, match, sPred, sAct

#main(useCOO=False, verbose=True)
