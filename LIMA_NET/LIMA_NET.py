import torch
import torch.nn as nn
import torch.nn.functional as F
import math


class LIMANET():
    def __init__(self, model, inputsize, dataloader):
        self.model = model

        self.dataloader = dataloader
        self.trainloader = dataloader.trainloader
        self.valloader = dataloader.valloader

        #self.loss = self.calcLoss()
        #self.loss = torch.nn.MSELoss()
        #self.loss = torch.nn.L1Loss()
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.0001)


        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(self.device)
        self.model = self.model.to(self.device)


    def train(self, n_epochs):
        print("Initial val loss: ", self.validate())
        for e in range(n_epochs):
            train_loss = self.train_one_epoch()
            val_loss = self.validate()
            print("Epoch ", e, " train-loss: ", train_loss, " val-loss: ", val_loss)


    def train_one_epoch(self):
        epoch_loss = 0
        for i, data in enumerate(self.trainloader):
            inputs, labels = data
            inputs = inputs.to(self.device)
            labels = labels.to(self.device)

            self.optimizer.zero_grad()

            outputs = self.model(inputs)
            #loss = self.loss(outputs, labels)
            loss = self.calcLoss(outputs, labels)
            loss.backward()

            self.optimizer.step()

            epoch_loss += loss.item()
        avg_loss = epoch_loss / (len(self.trainloader))

        return avg_loss


    def validate(self):
        self.model.train(False)
        vloss_total = 0


        for i, vdata in enumerate(self.valloader):
            vinputs, vlabels = vdata
            vinputs = vinputs.to(self.device)
            vlabels = vlabels.to(self.device)

            preds = self.model(vinputs)
            #vloss = self.loss(preds, vlabels)
            vloss = self.calcLoss(preds, vlabels)
            vloss_total += vloss.item()

            self.calcAccuracy(preds, vlabels)
        print(vloss_total)
        vloss_avg = vloss_total / (len(self.valloader) * 32)
        self.model.train(True)

        return vloss_avg

    def calcLoss(self, predictions, labels):
        error = predictions.sub(labels)
        error = torch.square(error)
        error = torch.sum(error, 1)
        #error = torch.sqrt(error)
        mean_err = torch.mean(error)
        return mean_err


    def calcAccuracy(self, predictions, labels):
        error = predictions.sub(labels)
        error = torch.square(error)
        MSE = torch.sum(error, 1)
        ME = torch.sqrt(error)

        min_err = torch.min(ME)
        max_err = torch.max(ME)
        mean_err = torch.mean(MSE)
        #print("Min ", min_err.item())
        print("Mean ", mean_err.item())
        #print("Max ", max_err.item())
        return min_err, max_err, mean_err