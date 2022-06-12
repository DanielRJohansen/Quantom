import torch
import torch.nn as nn
import torch.nn.functional as F

def makeBlock(inputs, outputs, size):
    return nn.Sequential(
        #nn.BatchNorm1d(inputs),
        nn.Linear(inputs, size),
        nn.ReLU(),

        nn.Linear(size, outputs),
    )

def makeBlock1(inputs, outputs, size):
    return nn.Sequential(
        nn.Linear(inputs, size),
        nn.BatchNorm1d(size),
        nn.Linear(size, size),
        nn.ReLU(),
        nn.Linear(size, size),
        nn.ReLU(),
    )




class LIMADNN4(nn.Module):
    def __init__(self, n_neighbors, n_bins, inputsize):
        super(LIMADNN4, self).__init__()

        #self.n_neighbors = n_neighbors

        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        self.out_bins = n_bins


        self.lstmblock = nn.LSTM(input_size=4, hidden_size=8, num_layers=2, batch_first=True)
        self.lstmout_fc = nn.Sequential(
            nn.Linear(8, 8),
            nn.ReLU(),
            nn.Linear(8, 8),
            nn.ReLU()
        )

        self.middleblock = nn.Sequential(
            nn.Linear(8,8),
            nn.ReLU(),
            nn.Linear(8, 8),
            nn.ReLU(),
        )

        self.meanpredict_block = nn.Sequential(
            nn.Linear(8,8),
            nn.ReLU(),
            nn.Linear(8, 8),
            nn.ReLU(),
            nn.Linear(8, self.out_bins),
            #nn.ReLU()
            nn.Sigmoid()
        )
        self.stdpredict_block = nn.Sequential(
            nn.Linear(8, 8),
            nn.ReLU(),
            nn.Linear(8, 8),
            nn.ReLU(),
            nn.Linear(8, self.out_bins),
            #nn.ReLU()
            nn.Sigmoid()
        )

        #self.mixer_fc = nn.Linear(12+6, 16)

    def calcElementLengths(self, x):
        squares = torch.square(x)
        sums = torch.sum(squares, dim=2)
        lengths = torch.sqrt(sums)
        return lengths

    def forward(self, x):
        p0data = x[:,:2,:]
        p0data = p0data.view(-1, 6)
        base = p0data[:,:3]


        #force_scalar = self.calcElementLengths(base[:, None, :])

        neighbordata = x[:,2:,:]
        nd_lengths = self.calcElementLengths(neighbordata)
        print(torch.mean(nd_lengths))
        print(torch.std(nd_lengths))
        exit()
        neighbordata = torch.cat((neighbordata, nd_lengths[:,:,None]), dim=2) # Trick to add dummy dimension, so we can concat later

        #x1 = self.p0block(p0data)
        #x1 = self.p0simple(force_scalar)

        _, (x2, _) = self.lstmblock(neighbordata)           # Extract only projection from final element in sequence
        x2 = x2[-1,:,:]                                     # Extract only proj from final LSTM in block
        x2 = self.lstmout_fc(x2)

        x2 = self.middleblock(x2)
        means = self.meanpredict_block(x2)
        stds = self.stdpredict_block(x2) + 0.01         # Std must never be 0

        out = torch.stack((means, stds), dim=2)


        #mixed = torch.cat((x1, x2), dim=1)

        #out = self.p0out(x1)
        #out = self.mixerblock(mixed)

#        out = F.softmax(out, dim=1)
        return out, base















