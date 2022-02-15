import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import numpy as nn

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
        print(model)
        self.total_params = sum(p.numel() for p in model.parameters())
        print("Total parameters: ", self.total_params)


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
        loss_total = 0
        acc_total = 0

        for i, data in enumerate(self.valloader):
            inputs, labels = data
            inputs = inputs.to(self.device)
            labels = labels.to(self.device)

            pred, base = self.model(inputs)

            loss = self.calcLoss(pred, base, labels)
            loss_total += loss.item()
            acc = self.calcAccuracy(pred, base, labels)
            acc_total += acc.item()



        loss_avg = (loss_total-loss.item()) / (len(self.valloader)-1)  # Remove the last batch, as it might not be full size
        acc_avg = (acc_total - acc.item()) / (len(self.valloader) - 1)  # Remove the last batch, as it might not be full size
        self.model.train(True)

        return loss_avg, acc_avg

    def calcEuclideanError(self, predictions, labels):
        error = predictions.sub(labels)
        # error = predictions.sub(true_dF)
        sq_e = torch.square(error)
        sum_sq_e = torch.sum(sq_e, 1)
        sum_e = torch.sqrt(sum_sq_e)
        return sum_e

    def calcLossScalar(self, base, labels):
        # Two ways of putting error in context of force
        true_dF = labels.sub(base)
        scalar = torch.sqrt(torch.sum(torch.square(true_dF)))  # Dunno which scalar to use, try both i guess

        # scalar = torch.sqrt(torch.sum(torch.square(base)))


        scalar = torch.max(scalar, torch.tensor(0.0001))  # Dont let scalar become too small, avoid exploding loss
        return scalar

    def calcLoss(self, predictions, base, labels):
        euclidean_errors = self.calcEuclideanError(predictions, labels)
        scalar = self.calcLossScalar(base, labels)
        sum_of_errors = euclidean_errors.div(scalar)

        mean_err = torch.mean(sum_of_errors)

        return mean_err


    def calcAccuracy(self, predictions, base, labels):
        euclidean_errors = self.calcEuclideanError(predictions, labels)
        scalar = self.calcLossScalar(base, labels)
        sum_of_errors = euclidean_errors.div(scalar)

        ones = torch.ones(sum_of_errors.shape, dtype=torch.float32).to(self.device)
        zeroes = torch.zeros(sum_of_errors.shape, dtype=torch.float32).to(self.device)
        unbound_acc = ones.sub(sum_of_errors)
        stack = torch.stack((unbound_acc, zeroes))
        bounded_acc = torch.max(stack, dim=0)[0]

        mean_acc = torch.mean(bounded_acc)
        return mean_acc