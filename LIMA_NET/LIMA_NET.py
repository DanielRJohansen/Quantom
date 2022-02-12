import torch
import torch.nn as nn
import torch.nn.functional as F



class LIMANET():
    def __init__(self, model, inputsize, dataloader):
        self.model = model

        self.dataloader = dataloader
        self.trainloader = dataloader.trainloader
        self.valloader = dataloader.valloader

        self.loss = torch.nn.MSELoss()
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001)

        self.n_epochs = 100

        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(self.device)
        self.model = self.model.to(self.device)


    def train(self):
        print("Initial val loss: ", self.validate())
        for e in range(self.n_epochs):
            train_loss = self.train_one_epoch()
            val_loss = self.validate()
            print("Epoch ", e, " train-loss: ", train_loss, " val-loss: ", val_loss)


    def train_one_epoch(self):
        epoch_loss = 0
        batch_cnt = 0
        for i, data in enumerate(self.trainloader):
            inputs, labels = data
            inputs = inputs.to(self.device)
            labels = labels.to(self.device)

            self.optimizer.zero_grad()

            outputs = self.model(inputs)
            loss = self.loss(outputs, labels)
            loss.backward()

            self.optimizer.step()

            epoch_loss += loss.item()
            batch_cnt += 1


        avg_loss = epoch_loss / (batch_cnt)

        return avg_loss

    def validate(self):
        self.model.train(False)
        vloss_total = 0
        batch_cnt = 0
        for i, vdata in enumerate(self.valloader):
            vinputs, vlabels = vdata
            vinputs = vinputs.to(self.device)
            vlabels = vlabels.to(self.device)

            voutputs = self.model(vinputs)
            vloss = self.loss(voutputs, vlabels)
            vloss_total += vloss.item()
            batch_cnt += 1

        vloss_avg = vloss_total / (batch_cnt)
        self.model.train(True)

        return vloss_avg

    #def calcAccuracy(self, ):