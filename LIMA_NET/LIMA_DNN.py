import torch
import torch.nn as nn
import torch.nn.functional as F


class LIMADNN(nn.Module):
    def __init__(self, inputsize, dataloader):
        super(LIMADNN, self).__init__()

        self.fc1 = nn.Linear(inputsize, inputsize)  # 5*5 from image dimension
        self.fc2 = nn.Linear(inputsize, 128)
        self.fc3 = nn.Linear(128, 16)
        self.fc4 = nn.Linear(16,3)

        self.dataloader = dataloader
        self.trainloader = dataloader.trainloader
        self.valloader = dataloader.valloader

        self.loss = torch.nn.MSELoss()
        self.optimizer = torch.optim.Adam(self.parameters(), lr=0.001)

        self.n_epochs = 100

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        x = self.fc4(x)
        return x
