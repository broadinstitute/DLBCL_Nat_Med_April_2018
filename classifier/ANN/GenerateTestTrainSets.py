import numpy as np
import pandas as pd
import os

def generateDFs(seed, originalReduced, usen1s = False, useProbs = False, COO=True):
    print(seed)
    np.random.seed(seed)
    if originalReduced:
        if not usen1s:
            reducedDF = pd.read_csv("DataTables/reducedDF.txt", sep = '\t')
        else:
            reducedDF = pd.read_csv("DataTables/reducedDFn1s.txt", sep = '\t')
    else:
        reducedDF = pd.read_csv("DataTables/fixedReducedCOO.txt", sep = '\t')
        if not COO:
            reducedDF.drop(reducedDF.columns[len(reducedDF.columns) - 1], axis=1, inplace=True)

    trainTest = pd.read_csv("DataTables/DLBCL_train_test_sets_01May2018.txt", sep = '\t')
    labels = pd.read_csv("DataTables/DLBCL_no_21q22.3.bestclus.txt", sep = '\t')
    probs = pd.read_csv("DataTables/heatmapMatrixAlpha.0.9.Power1.2.constantScalar.txt", sep = '\t')


    trainTest = trainTest[trainTest['pair_id'].isin(reducedDF.index)]
    trainSet = trainTest.loc[trainTest['train/test'] == "train"]
    trainSet = trainSet.reset_index(drop = True)
    trainSet = trainSet.sample(frac=1, random_state=seed)

    c1Subset = labels[labels['cluster'] == 1]
    c2Subset = labels[labels['cluster'] == 2]
    c3Subset = labels[labels['cluster'] == 3]
    c4Subset = labels[labels['cluster'] == 4]
    c5Subset = labels[labels['cluster'] == 5]

    c1Train = trainSet[trainSet['pair_id'].isin(c1Subset['SampleName'])]
    c2Train = trainSet[trainSet['pair_id'].isin(c2Subset['SampleName'])]
    c3Train = trainSet[trainSet['pair_id'].isin(c3Subset['SampleName'])]
    c4Train = trainSet[trainSet['pair_id'].isin(c4Subset['SampleName'])]
    c5Train = trainSet[trainSet['pair_id'].isin(c5Subset['SampleName'])]

    c1Train = c1Train.reindex(np.random.permutation(c1Train.index))
    c1Train = c1Train.reset_index(drop = True)
    c2Train = c2Train.reindex(np.random.permutation(c2Train.index))
    c2Train = c2Train.reset_index(drop=True)
    c3Train = c3Train.reindex(np.random.permutation(c3Train.index))
    c3Train = c3Train.reset_index(drop=True)
    c4Train = c4Train.reindex(np.random.permutation(c4Train.index))
    c4Train = c4Train.reset_index(drop=True)
    c5Train = c5Train.reindex(np.random.permutation(c5Train.index))
    c5Train = c5Train.reset_index(drop=True)

    validationSize = 0.25
    cutidx = int(c1Train.shape[0] * validationSize)
    validationC1 = c1Train.iloc[0:cutidx, :]
    cutidx = int(c2Train.shape[0] * validationSize)
    validationC2 = c2Train.iloc[0:cutidx, :]
    cutidx = int(c3Train.shape[0] * validationSize)
    validationC3 = c3Train.iloc[0:cutidx, :]
    cutidx = int(c4Train.shape[0] * validationSize)
    validationC4 = c4Train.iloc[0:cutidx, :]
    cutidx = int(c5Train.shape[0] * validationSize)
    validationC5 = c5Train.iloc[0:cutidx, :]

    fullValidationSet = validationC1.append(validationC2).append(validationC3).append(validationC4).append(validationC5)
    fullValidationSet = fullValidationSet.reset_index(drop=True)

    fullTrainSet = trainSet[~trainSet['pair_id'].isin(fullValidationSet['pair_id'])]

    validationFrame = reducedDF[reducedDF.index.isin(fullValidationSet['pair_id'])]
    trainFrame = reducedDF[reducedDF.index.isin(fullTrainSet['pair_id'])]

    #trainFrame = trainFrame.sample(frac=1, random_state=seed)

    appendTrainLabels = []
    appendValidationLabels = []
    for i in range(0,trainFrame.shape[0]):
        sampleName = trainFrame.iloc[i, :].name
        if not useProbs:
            assignedCluster = labels.loc[labels['SampleName'] == sampleName,'cluster']
            appendTrainLabels.append(assignedCluster.values[0])
        else:
            currProbs = probs.loc[probs['sample'] == sampleName, ['P1', 'P2', 'P3', 'P4', 'P5']]
            appVal = list(currProbs.values[0])
            appendTrainLabels.append(appVal)

    for i in range(0, validationFrame.shape[0]):
        sampleName = validationFrame.iloc[i,:].name
        if not useProbs:
            assignedCluster = labels.loc[labels['SampleName'] == sampleName, 'cluster']
            appendValidationLabels.append(assignedCluster.values[0])
        else:
            currProbs = probs.loc[probs['sample'] == sampleName,['P1','P2','P3','P4','P5']]
            appVal = list(currProbs.values[0])
            appendValidationLabels.append(appVal)

    return trainFrame,validationFrame,appendTrainLabels,appendValidationLabels

def generateExpandedSets(seed = np.random.randint(low=0, high = 99999999), expandSize=100):
    np.random.seed(seed)
    reducedDF = pd.read_csv("DataTables/fixedReducedCOO.txt", sep='\t')
    trainTest = pd.read_csv("DataTables/DLBCL_train_test_sets_01May2018.txt", sep='\t')
    labels = pd.read_csv("DataTables/DLBCL_no_21q22.3.bestclus.txt", sep='\t')
    probs = pd.read_csv("DataTables/baselineProbTable.txt", sep='\t', header=0)
    trainTest = trainTest[trainTest['pair_id'].isin(reducedDF.index)]
    trainSet = trainTest.loc[trainTest['train/test'] == "train"]
    trainSet = trainSet.reset_index(drop=True)

    c1Subset = labels[labels['cluster'] == 1]
    c2Subset = labels[labels['cluster'] == 2]
    c3Subset = labels[labels['cluster'] == 3]
    c4Subset = labels[labels['cluster'] == 4]
    c5Subset = labels[labels['cluster'] == 5]

    c1Train = trainSet[trainSet['pair_id'].isin(c1Subset['SampleName'])]
    c2Train = trainSet[trainSet['pair_id'].isin(c2Subset['SampleName'])]
    c3Train = trainSet[trainSet['pair_id'].isin(c3Subset['SampleName'])]
    c4Train = trainSet[trainSet['pair_id'].isin(c4Subset['SampleName'])]
    c5Train = trainSet[trainSet['pair_id'].isin(c5Subset['SampleName'])]

    c1Train = c1Train.reindex(np.random.permutation(c1Train.index))
    c1Train = c1Train.reset_index(drop=True)
    c2Train = c2Train.reindex(np.random.permutation(c2Train.index))
    c2Train = c2Train.reset_index(drop=True)
    c3Train = c3Train.reindex(np.random.permutation(c3Train.index))
    c3Train = c3Train.reset_index(drop=True)
    c4Train = c4Train.reindex(np.random.permutation(c4Train.index))
    c4Train = c4Train.reset_index(drop=True)
    c5Train = c5Train.reindex(np.random.permutation(c5Train.index))
    c5Train = c5Train.reset_index(drop=True)

    validationSize = .20
    cutidx = int(c1Train.shape[0] * validationSize)
    validationC1 = c1Train.iloc[0:cutidx, :]
    cutidx = int(c2Train.shape[0] * validationSize)
    validationC2 = c2Train.iloc[0:cutidx, :]
    cutidx = int(c3Train.shape[0] * validationSize)
    validationC3 = c3Train.iloc[0:cutidx, :]
    cutidx = int(c4Train.shape[0] * validationSize)
    validationC4 = c4Train.iloc[0:cutidx, :]
    cutidx = int(c5Train.shape[0] * validationSize)
    validationC5 = c5Train.iloc[0:cutidx, :]

    fullValidationSet = validationC1.append(validationC2).append(validationC3).append(validationC4).append(validationC5)
    fullValidationSet = fullValidationSet.reset_index(drop=True)

    fullTrainSet = trainSet[~trainSet['pair_id'].isin(fullValidationSet['pair_id'])]

    validationFrame = reducedDF[reducedDF.index.isin(fullValidationSet['pair_id'])]
    trainFrame = reducedDF[reducedDF.index.isin(fullTrainSet['pair_id'])]

    expandedDF = trainFrame
    expandedDF = expandedDF.append([expandedDF]*20)

    expandedDF = expandedDF.sample(frac=1, random_state=seed)
    shuffledCols = expandedDF.columns.values
    np.random.shuffle(shuffledCols)
    expandedDF = expandedDF[trainFrame.columns]
    validationFrame = validationFrame[trainFrame.columns]


    appendTrainLabels = []
    appendValidationLabels = []
    for i in range(0, expandedDF.shape[0]):
        sampleName = expandedDF.iloc[i, :].name
        currProb = probs.loc[sampleName]
        currProb = np.asarray(currProb.iloc[0:5])
        drawnCluster = np.random.choice([1, 2, 3, 4, 5], p=currProb)
        appendTrainLabels.append(drawnCluster)

        # sampleName = expandedDF.iloc[i, :].name
        # assignedCluster = labels.loc[labels['SampleName'] == sampleName, 'cluster']
        # appendTrainLabels.append(assignedCluster.values[0])

    for i in range(0, validationFrame.shape[0]):
        sampleName = validationFrame.iloc[i, :].name
        assignedCluster = labels.loc[labels['SampleName'] == sampleName, 'cluster']
        appendValidationLabels.append(assignedCluster.values[0])

    return expandedDF, validationFrame, appendTrainLabels, appendValidationLabels


def generateLowExpandedDF(seed, originalReduced, numExpand, useProbs = False, COO=True):
    print(seed)
    np.random.seed(seed)
    if originalReduced:
        reducedDF = pd.read_csv("DataTables/reducedDF.txt", sep = '\t')
    else:
        reducedDF = pd.read_csv("DataTables/fixedReducedCOO.txt", sep = '\t')
        if not COO:
            reducedDF.drop(reducedDF.columns[len(reducedDF.columns) - 1], axis=1, inplace=True)

    trainTest = pd.read_csv("DataTables/DLBCL_train_test_sets_01May2018.txt", sep = '\t')
    labels = pd.read_csv("DataTables/DLBCL_no_21q22.3.bestclus.txt", sep = '\t')
    probs = pd.read_csv("DataTables/heatmapMatrixAlpha.0.9.Power1.2.constantScalar.txt", sep = '\t')


    trainTest = trainTest[trainTest['pair_id'].isin(reducedDF.index)]
    trainSet = trainTest.loc[trainTest['train/test'] == "train"]
    trainSet = trainSet.reset_index(drop = True)

    c1Subset = labels[labels['cluster'] == 1]
    c2Subset = labels[labels['cluster'] == 2]
    c3Subset = labels[labels['cluster'] == 3]
    c4Subset = labels[labels['cluster'] == 4]
    c5Subset = labels[labels['cluster'] == 5]

    c1Train = trainSet[trainSet['pair_id'].isin(c1Subset['SampleName'])]
    c2Train = trainSet[trainSet['pair_id'].isin(c2Subset['SampleName'])]
    c3Train = trainSet[trainSet['pair_id'].isin(c3Subset['SampleName'])]
    c4Train = trainSet[trainSet['pair_id'].isin(c4Subset['SampleName'])]
    c5Train = trainSet[trainSet['pair_id'].isin(c5Subset['SampleName'])]

    c1Train = c1Train.reindex(np.random.permutation(c1Train.index))
    c1Train = c1Train.reset_index(drop = True)
    c2Train = c2Train.reindex(np.random.permutation(c2Train.index))
    c2Train = c2Train.reset_index(drop=True)
    c3Train = c3Train.reindex(np.random.permutation(c3Train.index))
    c3Train = c3Train.reset_index(drop=True)
    c4Train = c4Train.reindex(np.random.permutation(c4Train.index))
    c4Train = c4Train.reset_index(drop=True)
    c5Train = c5Train.reindex(np.random.permutation(c5Train.index))
    c5Train = c5Train.reset_index(drop=True)

    validationSize = .20
    cutidx = int(c1Train.shape[0] * validationSize)
    validationC1 = c1Train.iloc[0:cutidx, :]
    cutidx = int(c2Train.shape[0] * validationSize)
    validationC2 = c2Train.iloc[0:cutidx, :]
    cutidx = int(c3Train.shape[0] * validationSize)
    validationC3 = c3Train.iloc[0:cutidx, :]
    cutidx = int(c4Train.shape[0] * validationSize)
    validationC4 = c4Train.iloc[0:cutidx, :]
    cutidx = int(c5Train.shape[0] * validationSize)
    validationC5 = c5Train.iloc[0:cutidx, :]

    fullValidationSet = validationC1.append(validationC2).append(validationC3).append(validationC4).append(validationC5)
    fullValidationSet = fullValidationSet.reset_index(drop=True)

    fullTrainSet = trainSet[~trainSet['pair_id'].isin(fullValidationSet['pair_id'])]

    validationFrame = reducedDF[reducedDF.index.isin(fullValidationSet['pair_id'])]
    trainFrame = reducedDF[reducedDF.index.isin(fullTrainSet['pair_id'])]

    trainProbs = probs.loc[probs['sample'].isin(trainFrame.index.values)]
    trainProbs = trainProbs.reset_index(drop=True)
    lowConfs = trainProbs['maxprob'] < .85
    lowConfs = trainProbs[lowConfs]
    vals = lowConfs['sample'].values
    lowSubset = trainFrame.filter(items = vals, axis=0)
    toAppend = [lowSubset]*numExpand
    if toAppend:
        trainFrame = trainFrame.append(toAppend)

    trainFrame = trainFrame.sample(frac=1, random_state=seed)



    appendTrainLabels = []
    appendValidationLabels = []
    for i in range(0,trainFrame.shape[0]):
        sampleName = trainFrame.iloc[i, :].name
        if not useProbs:
            assignedCluster = labels.loc[labels['SampleName'] == sampleName,'cluster']
            appendTrainLabels.append(assignedCluster.values[0])
        else:
            currProbs = probs.loc[probs['sample'] == sampleName, ['P1', 'P2', 'P3', 'P4', 'P5']]
            appVal = list(currProbs.values[0])
            appendTrainLabels.append(appVal)

    for i in range(0, validationFrame.shape[0]):
        sampleName = validationFrame.iloc[i,:].name
        if not useProbs:
            assignedCluster = labels.loc[labels['SampleName'] == sampleName, 'cluster']
            appendValidationLabels.append(assignedCluster.values[0])
        else:
            currProbs = probs.loc[probs['sample'] == sampleName,['P1','P2','P3','P4','P5']]
            appVal = list(currProbs.values[0])
            appendValidationLabels.append(appVal)

    return trainFrame,validationFrame,appendTrainLabels,appendValidationLabels

def generateHighExpandedDF(seed, originalReduced, numExpand, useProbs = False, COO=True):
    print(seed)
    np.random.seed(seed)
    if originalReduced:
        reducedDF = pd.read_csv("DataTables/reducedDF.txt", sep = '\t')
    else:
        reducedDF = pd.read_csv("DataTables/fixedReducedCOO.txt", sep = '\t')
        if not COO:
            reducedDF.drop(reducedDF.columns[len(reducedDF.columns) - 1], axis=1, inplace=True)

    trainTest = pd.read_csv("DataTables/DLBCL_train_test_sets_01May2018.txt", sep = '\t')
    labels = pd.read_csv("DataTables/DLBCL_no_21q22.3.bestclus.txt", sep = '\t')
    probs = pd.read_csv("DataTables/heatmapMatrixAlpha.0.9.Power1.2.constantScalar.txt", sep = '\t')


    trainTest = trainTest[trainTest['pair_id'].isin(reducedDF.index)]
    trainSet = trainTest.loc[trainTest['train/test'] == "train"]
    trainSet = trainSet.reset_index(drop = True)

    c1Subset = labels[labels['cluster'] == 1]
    c2Subset = labels[labels['cluster'] == 2]
    c3Subset = labels[labels['cluster'] == 3]
    c4Subset = labels[labels['cluster'] == 4]
    c5Subset = labels[labels['cluster'] == 5]

    c1Train = trainSet[trainSet['pair_id'].isin(c1Subset['SampleName'])]
    c2Train = trainSet[trainSet['pair_id'].isin(c2Subset['SampleName'])]
    c3Train = trainSet[trainSet['pair_id'].isin(c3Subset['SampleName'])]
    c4Train = trainSet[trainSet['pair_id'].isin(c4Subset['SampleName'])]
    c5Train = trainSet[trainSet['pair_id'].isin(c5Subset['SampleName'])]

    c1Train = c1Train.reindex(np.random.permutation(c1Train.index))
    c1Train = c1Train.reset_index(drop = True)
    c2Train = c2Train.reindex(np.random.permutation(c2Train.index))
    c2Train = c2Train.reset_index(drop=True)
    c3Train = c3Train.reindex(np.random.permutation(c3Train.index))
    c3Train = c3Train.reset_index(drop=True)
    c4Train = c4Train.reindex(np.random.permutation(c4Train.index))
    c4Train = c4Train.reset_index(drop=True)
    c5Train = c5Train.reindex(np.random.permutation(c5Train.index))
    c5Train = c5Train.reset_index(drop=True)

    validationSize = .2
    cutidx = int(c1Train.shape[0] * validationSize)
    validationC1 = c1Train.iloc[0:cutidx, :]
    cutidx = int(c2Train.shape[0] * validationSize)
    validationC2 = c2Train.iloc[0:cutidx, :]
    cutidx = int(c3Train.shape[0] * validationSize)
    validationC3 = c3Train.iloc[0:cutidx, :]
    cutidx = int(c4Train.shape[0] * validationSize)
    validationC4 = c4Train.iloc[0:cutidx, :]
    cutidx = int(c5Train.shape[0] * validationSize)
    validationC5 = c5Train.iloc[0:cutidx, :]

    fullValidationSet = validationC1.append(validationC2).append(validationC3).append(validationC4).append(validationC5)
    fullValidationSet = fullValidationSet.reset_index(drop=True)

    fullTrainSet = trainSet[~trainSet['pair_id'].isin(fullValidationSet['pair_id'])]

    validationFrame = reducedDF[reducedDF.index.isin(fullValidationSet['pair_id'])]
    trainFrame = reducedDF[reducedDF.index.isin(fullTrainSet['pair_id'])]


    trainProbs = probs.loc[probs['sample'].isin(trainFrame.index.values)]
    trainProbs = trainProbs.reset_index(drop=True)
    highConfs = trainProbs['maxprob'] > .94
    highConfs = trainProbs[highConfs]
    vals = highConfs['sample'].values
    lowSubset = trainFrame.filter(items = vals, axis=0)
    toAppend = [lowSubset]*numExpand
    if toAppend:
        trainFrame = trainFrame.append(toAppend)

    trainFrame = trainFrame.sample(frac=1, random_state=seed)

    appendTrainLabels = []
    appendValidationLabels = []
    for i in range(0,trainFrame.shape[0]):
        sampleName = trainFrame.iloc[i, :].name
        if not useProbs:
            assignedCluster = labels.loc[labels['SampleName'] == sampleName,'cluster']
            appendTrainLabels.append(assignedCluster.values[0])
        else:
            currProbs = probs.loc[probs['sample'] == sampleName, ['P1', 'P2', 'P3', 'P4', 'P5']]
            appVal = list(currProbs.values[0])
            appendTrainLabels.append(appVal)

    for i in range(0, validationFrame.shape[0]):
        sampleName = validationFrame.iloc[i,:].name
        if not useProbs:
            assignedCluster = labels.loc[labels['SampleName'] == sampleName, 'cluster']
            appendValidationLabels.append(assignedCluster.values[0])
        else:
            currProbs = probs.loc[probs['sample'] == sampleName,['P1','P2','P3','P4','P5']]
            appVal = list(currProbs.values[0])
            appendValidationLabels.append(appVal)

    return trainFrame,validationFrame,appendTrainLabels,appendValidationLabels

def generateFullDF(seed, useProbs = False):

    np.random.seed(seed)
    fullDF = pd.read_csv("DataTables/fullDF.txt", sep='\t')
    trainTest = pd.read_csv("DataTables/DLBCL_train_test_sets_01May2018.txt", sep='\t')
    labels = pd.read_csv("DataTables/DLBCL_no_21q22.3.bestclus.txt", sep='\t')
    probs = pd.read_csv("DataTables/heatmapMatrixAlpha.0.9.Power1.2.constantScalar.txt", sep='\t')

    trainTest = trainTest[trainTest['pair_id'].isin(fullDF.index)]
    trainSet = trainTest.loc[trainTest['train/test'] == "train"]
    trainSet = trainSet.reset_index(drop=True)

    c1Subset = labels[labels['cluster'] == 1]
    c2Subset = labels[labels['cluster'] == 2]
    c3Subset = labels[labels['cluster'] == 3]
    c4Subset = labels[labels['cluster'] == 4]
    c5Subset = labels[labels['cluster'] == 5]

    c1Train = trainSet[trainSet['pair_id'].isin(c1Subset['SampleName'])]
    c2Train = trainSet[trainSet['pair_id'].isin(c2Subset['SampleName'])]
    c3Train = trainSet[trainSet['pair_id'].isin(c3Subset['SampleName'])]
    c4Train = trainSet[trainSet['pair_id'].isin(c4Subset['SampleName'])]
    c5Train = trainSet[trainSet['pair_id'].isin(c5Subset['SampleName'])]

    c1Train = c1Train.reindex(np.random.permutation(c1Train.index))
    c1Train = c1Train.reset_index(drop=True)
    c2Train = c2Train.reindex(np.random.permutation(c2Train.index))
    c2Train = c2Train.reset_index(drop=True)
    c3Train = c3Train.reindex(np.random.permutation(c3Train.index))
    c3Train = c3Train.reset_index(drop=True)
    c4Train = c4Train.reindex(np.random.permutation(c4Train.index))
    c4Train = c4Train.reset_index(drop=True)
    c5Train = c5Train.reindex(np.random.permutation(c5Train.index))
    c5Train = c5Train.reset_index(drop=True)

    validationSize = .20
    cutidx = int(c1Train.shape[0] * validationSize)
    validationC1 = c1Train.iloc[0:cutidx, :]
    cutidx = int(c2Train.shape[0] * validationSize)
    validationC2 = c2Train.iloc[0:cutidx, :]
    cutidx = int(c3Train.shape[0] * validationSize)
    validationC3 = c3Train.iloc[0:cutidx, :]
    cutidx = int(c4Train.shape[0] * validationSize)
    validationC4 = c4Train.iloc[0:cutidx, :]
    cutidx = int(c5Train.shape[0] * validationSize)
    validationC5 = c5Train.iloc[0:cutidx, :]

    fullValidationSet = validationC1.append(validationC2).append(validationC3).append(validationC4).append(validationC5)
    fullValidationSet = fullValidationSet.reset_index(drop=True)

    fullTrainSet = trainSet[~trainSet['pair_id'].isin(fullValidationSet['pair_id'])]

    validationFrame = fullDF[fullDF.index.isin(fullValidationSet['pair_id'])]
    trainFrame = fullDF[fullDF.index.isin(fullTrainSet['pair_id'])]

    trainFrame = trainFrame.sample(frac=1, random_state=seed)

    appendTrainLabels = []
    appendValidationLabels = []
    for i in range(0, trainFrame.shape[0]):
        sampleName = trainFrame.iloc[i, :].name
        if not useProbs:
            assignedCluster = labels.loc[labels['SampleName'] == sampleName, 'cluster']
            appendTrainLabels.append(assignedCluster.values[0])
        else:
            currProbs = probs.loc[probs['sample'] == sampleName, ['P1', 'P2', 'P3', 'P4', 'P5']]
            appVal = list(currProbs.values[0])
            appendTrainLabels.append(appVal)

    for i in range(0, validationFrame.shape[0]):
        sampleName = validationFrame.iloc[i, :].name
        if not useProbs:
            assignedCluster = labels.loc[labels['SampleName'] == sampleName, 'cluster']
            appendValidationLabels.append(assignedCluster.values[0])
        else:
            currProbs = probs.loc[probs['sample'] == sampleName, ['P1', 'P2', 'P3', 'P4', 'P5']]
            appVal = list(currProbs.values[0])
            appendValidationLabels.append(appVal)

    return trainFrame, validationFrame, appendTrainLabels, appendValidationLabels


def generateKFoldDFs(seed, originalReduced, folds, COO, downSample, useSV,
                     useCNA, GD, useProbs, fixedValidationSets, fullFeatures):
    print(seed)
    print(os.getcwd())
    np.random.seed(seed)
    if originalReduced:
        if not useSV and not useCNA:
            reducedDF = pd.read_csv("DataTables/reducedDF_noSV_noCNA.tsv", sep='\t')
        elif not useCNA:
            reducedDF = pd.read_csv("DataTables/reducedDF_noCNA.tsv", sep='\t')
        elif not useSV:
            reducedDF = pd.read_csv("DataTables/reducedDF_noSV.tsv", sep='\t')
        elif GD:
            reducedDF = pd.read_csv("DataTables/reducedDF_GD.txt", sep='\t')
        else:
            reducedDF = pd.read_csv("DataTables/reducedDF.txt", sep = '\t')
        if downSample:
            fixed = pd.read_csv("DataTables/fixedReducedCOO.txt", sep = '\t')
            reducedDF = reducedDF[reducedDF.index.isin(fixed.index)]
    else:
        reducedDF = pd.read_csv("DataTables/fixedReducedCOO.txt", sep = '\t')
        if not COO:
            reducedDF.drop(reducedDF.columns[len(reducedDF.columns) - 1], axis=1, inplace=True)

    if fullFeatures:
        reducedDF = pd.read_csv("DataTables/FullDF161.tsv", sep='\t')

    print(reducedDF.shape)
    trainTest = pd.read_csv("DataTables/DLBCL_train_test_sets_01May2018.txt", sep = '\t')
    labels = pd.read_csv("DataTables/DLBCL_no_21q22.3.bestclus.txt", sep = '\t')
    if useProbs:
        probs = pd.read_csv("DataTables/heatmapMatrixAlpha.0.9.Power1.2.constantScalar.txt", sep = '\t')


    trainTest = trainTest[trainTest['pair_id'].isin(reducedDF.index)]
    trainSet = trainTest.loc[trainTest['train/test'] == "train"]
    print(trainSet.shape)
    trainSet = trainSet.reset_index(drop = True)
    trainSet = trainSet.sample(frac=1, random_state=seed)

    c1Subset = labels[labels['cluster'] == 1]
    c2Subset = labels[labels['cluster'] == 2]
    c3Subset = labels[labels['cluster'] == 3]
    c4Subset = labels[labels['cluster'] == 4]
    c5Subset = labels[labels['cluster'] == 5]

    c1Train = trainSet[trainSet['pair_id'].isin(c1Subset['SampleName'])]
    c2Train = trainSet[trainSet['pair_id'].isin(c2Subset['SampleName'])]
    c3Train = trainSet[trainSet['pair_id'].isin(c3Subset['SampleName'])]
    c4Train = trainSet[trainSet['pair_id'].isin(c4Subset['SampleName'])]
    c5Train = trainSet[trainSet['pair_id'].isin(c5Subset['SampleName'])]

    c1Train = c1Train.reindex(np.random.permutation(c1Train.index))
    c1Train = c1Train.reset_index(drop = True)
    c2Train = c2Train.reindex(np.random.permutation(c2Train.index))
    c2Train = c2Train.reset_index(drop=True)
    c3Train = c3Train.reindex(np.random.permutation(c3Train.index))
    c3Train = c3Train.reset_index(drop=True)
    c4Train = c4Train.reindex(np.random.permutation(c4Train.index))
    c4Train = c4Train.reset_index(drop=True)
    c5Train = c5Train.reindex(np.random.permutation(c5Train.index))
    c5Train = c5Train.reset_index(drop=True)

    c1Samples = c1Train['pair_id'].tolist()
    c2Samples = c2Train['pair_id'].tolist()
    c3Samples = c3Train['pair_id'].tolist()
    c4Samples = c4Train['pair_id'].tolist()
    c5Samples = c5Train['pair_id'].tolist()

    validationSets = [[] for k in range(folds)]
    for i in range(0, len(c1Samples)):
        idx = i % 4
        validationSets[idx].append(c1Samples[i])
    for i in range(0, len(c2Samples)):
        idx = i % 4
        validationSets[idx].append(c2Samples[i])
    for i in range(0, len(c3Samples)):
        idx = i % 4
        validationSets[idx].append(c3Samples[i])
    for i in range(0, len(c4Samples)):
        idx = i % 4
        validationSets[idx].append(c4Samples[i])
    for i in range(0, len(c5Samples)):
        idx = i % 4
        validationSets[idx].append(c5Samples[i])

    #overwrite the validation frames if fixed frames. Not efficient. Not the philosophically correct way. Sorry.
    #hard coded in 4 folds too.

    if fixedValidationSets:

        v1 = pd.read_csv("DataTables/validationSet.1.txt", sep='\t', header=None)
        v2 = pd.read_csv("DataTables/validationSet.2.txt", sep='\t', header=None)
        v3 = pd.read_csv("DataTables/validationSet.3.txt", sep='\t', header=None)
        v4 = pd.read_csv("DataTables/validationSet.4.txt", sep='\t', header=None)

        validationSets = [v1.iloc[:,0].values, v2.iloc[:,0].values,v3.iloc[:,0].values,v4.iloc[:,0].values]



    trainingSets = []
    for k in range(0, folds):
        currVal = validationSets[k]
        tmp = trainSet[~trainSet['pair_id'].isin(currVal)]
        tmp = tmp.reset_index(drop=True)
        trainingSets.append(tmp['pair_id'].tolist())


    trainFrames = []
    validationFrames = []
    trainLabels = [[] for k in range(folds)]
    validationLabels = [[] for k in range(folds)]
    for k in range(0, folds):
        toAppend = reducedDF[reducedDF.index.isin(trainingSets[k])]
        toAppend = toAppend.reindex(trainingSets[k])
        trainFrames.append(toAppend)
        toAppend2 = reducedDF[reducedDF.index.isin(validationSets[k])]
        toAppend2 = toAppend2.reindex(validationSets[k])
        validationFrames.append(toAppend2)
        for i in range(0, len(trainingSets[k])):
            sampleName = trainingSets[k][i]
            if useProbs:
                currProbs = probs.loc[probs['sample'] == sampleName, ['P1', 'P2', 'P3', 'P4', 'P5']]
                appVal = list(currProbs.values[0])
            else:
                currLabel = labels.loc[labels['SampleName'] == sampleName, "cluster"]
                appVal = list(currLabel-1)
            trainLabels[k].append(appVal)
        for i in range(0, len(validationSets[k])):
            sampleName = validationSets[k][i]
            if useProbs:
                currProbs = probs.loc[probs['sample'] == sampleName, ['P1', 'P2', 'P3', 'P4', 'P5']]
                appVal = list(currProbs.values[0])
            else:
                currLabel = labels.loc[labels['SampleName'] == sampleName, "cluster"]
                appVal = list(currLabel - 1)
            validationLabels[k].append(appVal)



    return trainFrames, validationFrames, trainLabels, validationLabels

#generateKFoldDFs(123, True, 4)
#generateDFs(123, True)
