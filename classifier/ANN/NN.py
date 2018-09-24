import torch
import torch.nn as nn
import torch.nn.functional as F

class Net(nn.Module):

    def __init__(self, hlSize,inputSize,outputSize):
        super(Net, self).__init__()
        self.inputLayer = nn.Linear(inputSize, hlSize)
        self.hl1 = nn.Linear(hlSize, hlSize)
        self.bias = nn.Parameter(torch.ones(1))
        self.lastBias = nn.Parameter(torch.ones(1))
        #self.hl2 = nn.Linear(hlSize, hlSize)
        self.outputLayer = nn.Linear(hlSize, outputSize)

    def forward(self, x):
        m = nn.Tanh()
        # x = F.relu(self.inputLayer(x))
        # x = F.relu(self.hl1(x))
        #x = F.relu(self.hl2(x))
        x = self.inputLayer(x)
        x = m(x)
        x = self.hl1(x) + self.bias
        x = m(x)
        x = self.outputLayer(x) + self.lastBias
        x = F.softmax(x)
        # m = nn.Sigmoid()
        # x = m(x)
        # x = x/torch.sum(x)
        return x