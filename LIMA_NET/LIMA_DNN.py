import torch
import torch.nn as nn
import torch.nn.functional as F


class LIMADNN(nn.Module):
    def __init__(self, n_neighbors):
        super(LIMADNN, self).__init__()

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

        self.fc_fullstack1 = nn.Linear(16, 8)

        self.fc_out = nn.Linear(8, 3)


    def __forward(self, x):
        force_prev = x[:,:3]
        npos = x[:,3:self.n_neighbors*3+3]
        nforce = x[:,self.n_neighbors*3+3:self.n_neighbors*3*2+3]

        npos = self.fc_npos1(npos)
        npos = self.fc_npos2(npos)
        npos = F.relu(npos)
        
        nforce = self.fc_nforce1(nforce)
        nforce = self.fc_nforce2(nforce)
        nforce = F.relu(nforce)
        
        posforce = torch.cat((npos, nforce), 1)
        posforce = self.fc_forcepos_stack_1(posforce)
        posforce = self.fc_forcepos_stack_2(posforce)

        sforce = self.fc_sforce1(force_prev)

        stack = torch.cat((posforce, sforce), 1)
        stack = self.fc_fullstack1(stack)

        out = self.fc_out(stack)
        out = torch.sigmoid(out)
        
        out = out.mul(force_prev)
        out = out.add(force_prev)
        return out, force_prev

    def __forward0neighbors(self, x):
        zeroes = torch.zeros((x.shape[0], 6), dtype=torch.float32).to(self.device)
        x = torch.cat((x, zeroes), dim=1)

        return self.__forward(x)
