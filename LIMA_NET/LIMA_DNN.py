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

    def train(self):
        for e in range(self.n_epochs):
            epoch_loss_avg = self.train_one_epoch()
            print("Epoch ", e, " loss: ", epoch_loss_avg)

    def train_one_epoch(self):
        epoch_loss = 0
        batch_cnt = 0
        # Here, we use enumerate(training_loader) instead of
        # iter(training_loader) so that we can track the batch
        # index and do some intra-epoch reporting
        for i, data in enumerate(self.trainloader):
            # Every data instance is an input + label pair
            inputs, labels = data

            # Zero your gradients for every batch!
            self.optimizer.zero_grad()

            # Make predictions for this batch
            outputs = self.forward(inputs)
            # Compute the loss and its gradients
            loss = self.loss(outputs, labels)
            loss.backward()

            # Adjust learning weights
            self.optimizer.step()

            # Gather data and report
            epoch_loss += loss.item()
            batch_cnt += 1
            if i % 10 == 9:
                #last_loss = running_loss / 1000  # loss per batch
                #print('  batch {} loss: {}'.format(i + 1, last_loss))
                #tb_x = epoch_index * len(self.training_loader) + i + 1
                #tb_writer.add_scalar('Loss/train', last_loss, tb_x)
                running_loss = 0.

        avg_loss = epoch_loss / (self.dataloader.n_train)

        return avg_loss