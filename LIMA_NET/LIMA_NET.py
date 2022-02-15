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

            preds, bases = self.model(inputs)
            #loss = self.loss(outputs, labels)
            loss = self.calcLoss(preds, bases, labels)
            loss.backward()

            self.optimizer.step()

            epoch_loss += loss.item()

        avg_loss = (epoch_loss-loss.item()) / (len(self.trainloader)-1) # Remove the last batch, as it might not be full size

        return avg_loss


    def validate(self):
        self.model.train(False)
        vloss_total = 0


        for i, vdata in enumerate(self.valloader):
            vinputs, vlabels = vdata
            vinputs = vinputs.to(self.device)
            labels = vlabels.to(self.device)

            pred, base = self.model(vinputs)
            #vloss = self.loss(preds, vlabels)
            vloss = self.calcLoss(pred, base, labels)
            vloss_total += vloss.item()

            self.calcAccuracy(pred, base, labels)

        #print("pred: ", pred[0])
        #print("label ", labels[0])

        vloss_avg = (vloss_total-vloss.item()) / (len(self.valloader)-1)  # Remove the last batch, as it might not be full size
        self.model.train(True)

        return vloss_avg

    def calcLoss(self, predictions, base, labels):
        true_dF = labels.sub(base)

        error = predictions.sub(labels)
        #error = predictions.sub(true_dF)
        sq_e = torch.square(error)
        sum_sq_e = torch.sum(sq_e, 1)
        sum_e = torch.sqrt(sum_sq_e)

        scalar = torch.sqrt(torch.sum(torch.square(true_dF)))  # Dunno which scalar to use, try both i guess
        # scalar = torch.sqrt(torch.sum(torch.square(base)))
        scalar = torch.max(scalar, torch.tensor(0.001))
        sum_e = sum_e.div(scalar)
        mean_err = torch.mean(sum_e)

        return mean_err


    def calcAccuracy(self, predictions, base, labels):
        true_dF = labels.sub(base)

        #error = predictions.sub(labels)
        error = predictions.sub(true_dF)
        sq_e = torch.square(error)
        sum_sq_e = torch.sum(sq_e, 1)
        sum_e = torch.sqrt(sum_sq_e)

        scalar = torch.sqrt(torch.sum(torch.square(true_dF)))   # Dunno which scalar to use, try both i guess
        #scalar = torch.sqrt(torch.sum(torch.square(base)))
        scalar = torch.min(scalar, torch.tensor(0.000001))
        sum_e = sum_e.div(scalar)


        min_err = torch.min(sum_e)
        max_err = torch.max(sum_e)
        mean_err = torch.mean(sum_e)
        #print("Min ", min_err.item())
        if (mean_err > 10000):
            pass
            #print("Mean ", mean_err.item())
            #print("E ", error)
            #print(predictions)
            #print(labels)

        #print("Max ", max_err.item())
        return min_err, max_err, mean_err