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
        #self.block = makeBlock((self.n_neighbors) * 3, self.out_bins, 64)

        #self.block = makeBlock1((self.n_neighbors) * 3, self.out_bins, 256)
        #self.block = makeBlock1(inputsize, self.out_bins, 64)
        self.p0block = makeBlock1(6, 8, 8)
        self.p0out = nn.Sequential(
            nn.Linear(32, self.out_bins),
            nn.ReLU()
        )

        self.p0simple = nn.Sequential(
            nn.Linear(1, 4),
            nn.ReLU(),
            nn.Linear(4, 8),
            nn.ReLU()
        )


        self.lstmblock = nn.LSTM(input_size=4, hidden_size=32, num_layers=4, batch_first=True, proj_size=12)
        self.lstmout_fc = nn.Sequential(
            nn.Linear(12, 8),
            nn.ReLU()
        )


        self.mixerblock = nn.Sequential(
            nn.Linear(8+8, 16),
            nn.ReLU(),
            nn.Linear(16, self.out_bins),
            nn.ReLU()
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
        force_scalar = self.calcElementLengths(base[:, None, :])

        neighbordata = x[:,2:,:]
        nd_lengths = self.calcElementLengths(neighbordata)
        neighbordata = torch.cat((neighbordata, nd_lengths[:,:,None]), dim=2) # Trick to add dummy dimension, so we can concat later



        #x1 = self.p0block(p0data)
        x1 = self.p0simple(force_scalar)

        _, (x2, _) = self.lstmblock(neighbordata)           # Extract only projection from final element in sequence
        x2 = x2[-1,:,:]                                     # Extract only proj from final LSTM in block
        x2 = self.lstmout_fc(x2)

        mixed = torch.cat((x1, x2), dim=1)

        #out = self.p0out(x1)
        out = self.mixerblock(mixed)

        out = F.softmax(out, dim=1)
        return out, base















class LIMADNN3(nn.Module):
    def __init__(self, n_neighbors):
        super(LIMADNN3, self).__init__()

        self.n_neighbors = n_neighbors
        self.forward = self.__forward



        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        self.state_size = 16
        self.block = makeBlock((self.n_neighbors+1)*12, 3, 32)



    def __forward(self, x):
        force_prev = x[:, 0, 6:9]
        forcechange_prev = x[:, 0, 9:12]

        x = x.view(x.shape[0], x.shape[1]*x.shape[2])

        change_pred = self.block(x)


        #change_pred = self.out_layer(x)
        out = force_prev.add(forcechange_prev.mul(change_pred))

        #out = self.out_layer(x)
        #out = out.add(force_prev)



        return out, force_prev+forcechange_prev










def makeAtomBlock(inputs, outputs):
    return nn.Sequential(
            nn.Linear(inputs, 16),
            #nn.BatchNorm1d(),
            nn.ReLU(),
            nn.Linear(16, 16),
            nn.ReLU(),
            nn.Linear(16, 32),
            nn.ReLU(),
            nn.Linear(32, 16),
            # nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.Linear(16, outputs),
            nn.ReLU()
        )

class LIMADNN2(nn.Module):
    def __init__(self, n_neighbors):
        super(LIMADNN2, self).__init__()

        if (n_neighbors > 0):  # Model can't handle 0 neighbors, so we fake and add a flaot3(0,0,0)
            self.n_neighbors = n_neighbors
            self.forward = self.__forward
        else:
            self.n_neighbors = 1
            self.forward = self.__forward0neighbors


        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        self.state_size = 16
        self.blockout_size = 8

#        self.query_atom_block = makeAtomBlock(12, self.state_size)

        self.atom_block = makeAtomBlock(12, self.blockout_size)

        self.state_layer = nn.Linear(self.blockout_size, self.state_size)
        self.out_layer = nn.Linear(self.state_size, 3)

    def __forward(self, x):
        print(x.shape)
        exit()
        force_prev = x[:, 0, 6:9]
        forcechange_prev = x[:, 0, 9:12]

        #state = self.query_atom_block(x[:,0,:])
        queryblock_out = self.atom_block(x[:,0,:])
        state = self.state_layer(queryblock_out)

        #state = torch.ones(self.state_size, dtype=torch.float32).to(self.device)

        for i in range(0, x.shape[1]):
            block_out = self.atom_block(x[ : , i, : ])
            state_addition = self.state_layer(block_out)
            state.add(state_addition)


        out = self.out_layer(state)

        return out, force_prev+forcechange_prev

    def __forward0neighbors(self, x):
        zeroes = torch.zeros((x.shape[0], 6), dtype=torch.float32).to(self.device)
        x = torch.cat((x, zeroes), dim=1)

        return self.__forward(x)


class LIMADNN1(nn.Module):
    def __init__(self, n_neighbors):
        super(LIMADNN1, self).__init__()

        if (n_neighbors > 0):              # Model can't handle 0 neighbors, so we fake and add a flaot3(0,0,0)
            self.n_neighbors = n_neighbors
            self.forward = self.__forward
        else:
            self.n_neighbors = 1
            self.forward = self.__forward0neighbors


        neighbors_out = min(self.n_neighbors * 2, 32)

        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


        self.fc_npos1 = nn.Linear(self.n_neighbors*3, self.n_neighbors*3)
        self.fc_npos2 = nn.Linear(self.n_neighbors * 3, neighbors_out)

        self.fc_nforce1 = nn.Linear(self.n_neighbors*3, self.n_neighbors*3)
        self.fc_nforce2 = nn.Linear(self.n_neighbors * 3, neighbors_out)

        self.fc_forcepos_stack_1 = nn.Linear(neighbors_out*2, neighbors_out*2)
        self.fc_forcepos_stack_2 = nn.Linear(neighbors_out*2, 16-3)

        self.fc_sforce1 = nn.Linear(3, 3)

        self.fc_fullstack1 = nn.Linear(16, 16)
        self.fc_fullstack2 = nn.Linear(16, 8)

        self.fc_out = nn.Linear(8, 3)


    def __forward(self, x):
        force_prev = x[:,:3]
        ones = torch.ones((x.shape[0],3)).to(self.device)

        npos = x[:,3:self.n_neighbors*3+3]
        nforce = x[:,self.n_neighbors*3+3:self.n_neighbors*3*2+3]

        npos = self.fc_npos1(npos)
        npos = F.relu(npos)
        npos = self.fc_npos2(npos)
        npos = F.relu(npos)
        
        nforce = self.fc_nforce1(nforce)
        nforce = F.relu(nforce)
        nforce = self.fc_nforce2(nforce)
        nforce = F.relu(nforce)
        
        posforce = torch.cat((npos, nforce), 1)
        posforce = self.fc_forcepos_stack_1(posforce)
        posforce = F.relu(posforce)
        posforce = self.fc_forcepos_stack_2(posforce)
        posforce = F.relu(posforce)

        #sforce = self.fc_sforce1(force_prev)
        sforce = self.fc_sforce1(ones)

        stack = torch.cat((posforce, sforce), 1)
        stack = self.fc_fullstack1(stack)
        stack = self.fc_fullstack2(stack)

        out = self.fc_out(stack)
        #out = torch.sigmoid(out)
        
        #out = out.mul(force_prev)

        #out = out.add(force_prev)
        return out, force_prev

    def __forward0neighbors(self, x):
        zeroes = torch.zeros((x.shape[0], 6), dtype=torch.float32).to(self.device)
        x = torch.cat((x, zeroes), dim=1)

        return self.__forward(x)



