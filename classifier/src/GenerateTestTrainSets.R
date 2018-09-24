trainTest = read.table("DataTables/DLBCL_train_test_sets_01May2018.txt", header = TRUE)
trainTest = data.frame(trainTest)
trainingSet = as.vector(trainTest[trainTest$train.test == "train","pair_id"])
testingSet = as.vector(trainTest[trainTest$train.test == "test","pair_id"])
