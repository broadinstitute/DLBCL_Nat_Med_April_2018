import NN
import GenerateTestTrainSets as gt
import numpy as np
import torch
import torch.nn as nn
import pandas as pd
import matplotlib.pyplot as pl
from scipy.stats.stats import pearsonr

def main(seed, folds, original=True, COO=False,
         verbose=False, downSample = False, useSV = True,
         useCNA = True, GD=False, useProbs=True,
         fixedValidationSets=False, fullFeatures=False):
    trainFrames, validationFrames, trainLabels, validationLabels = gt.generateKFoldDFs(seed,
                                                                                       original,
                                                                                       folds,
                                                                                       COO,
                                                                                       downSample,
                                                                                       useSV,
                                                                                       useCNA,
                                                                                       GD,
                                                                                       useProbs=useProbs,
                                                                                       fixedValidationSets=fixedValidationSets,
                                                                                       fullFeatures=fullFeatures)



    #torch.manual_seed(seed)

    nets = []
    for k in range(0, folds):
        net = NN.Net(10, trainFrames[k].shape[1], 5)
        nets.append(net)
        print("Training size: " + str(trainFrames[k].shape[0]))

    trainValues = []
    validationValues = []
    sampleNames = []
    for k in range(0, folds):
        trainValues.append(trainFrames[k].values)
        validationValues.append(validationFrames[k].values)
        sampleNames.append(trainFrames[k].index.values)


    # lr = 0.01
    # num_epoch = 200

    if useProbs:
        num_epoch = 30
        lr = 0.1
    else:
        num_epoch = 100
        lr = 0.01

    if useProbs:
        MSE = True
    else:
        MSE = False

    if MSE:
        criterion = nn.MSELoss()
    else:
        criterion = nn.CrossEntropyLoss()

    optimizers = []
    for k in range(0, folds):
        optimizer = torch.optim.SGD(nets[k].parameters(), lr=lr, weight_decay=.001)
        optimizers.append(optimizer)

    trainInputsArr = [[] for k in range(folds)]
    trainLabelsArr = [[] for k in range(folds)]

    for k in range(0, folds):
        for i in range(0,len(trainValues[k])):
            if not MSE:
                toAppend = torch.tensor([trainValues[k][i]], requires_grad=True).float()
                toAppend2 = torch.tensor(trainLabels[k][i]).long()
            else:
                toAppend = torch.tensor([trainValues[k][i]], requires_grad=True).float()
                toAppend2 = torch.tensor([trainLabels[k][i]]).float()
            trainLabelsArr[k].append(toAppend2)
            trainInputsArr[k].append(toAppend)

    finalPredictions = []
    finalTruths = []
    trainedNets = []
    for k in range(0, folds):
        net = nets[k]
        optimizer = optimizers[k]
        currTrainInputs = trainInputsArr[k]
        currTrainLabels = trainLabelsArr[k]
        for epoch in range(0, num_epoch):
            avgloss = 0
            currPreds = []
            currTruths = []
            for i in range(0, len(currTrainInputs)):
                X = currTrainInputs[i]
                Y = currTrainLabels[i]
                #feedforward - backprop
                optimizer.zero_grad()
                out = net.forward(X)
                loss = criterion(out, Y)
                avgloss += loss.item()
                loss.backward()
                optimizer.step()
                if MSE:
                    currp = torch.max(Y, 1)[0].item()
                else:
                    currp = Y.item()
                #idx = torch.max(Y, 1)[1].item()
                if MSE:
                    val = np.max(out.data[0].numpy())
                else:
                    val = np.argmax(out.data[0].numpy())
                #val = val[idx]
                currPreds.append(val)
                currTruths.append(currp)

            avgloss = avgloss/len(trainInputsArr[k])
            finalPredictions = currPreds
            finalTruths = currTruths

            toprint = True
            if(toprint):
                if (epoch) % (num_epoch/10) == 0:
                    print ('Epoch [%d/%d] Loss: %.4f' % (epoch, num_epoch, avgloss))
        trainedNets.append(net)


    return nets, validationFrames, validationLabels