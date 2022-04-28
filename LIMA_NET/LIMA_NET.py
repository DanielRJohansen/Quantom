import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import numpy as np

class LIMANET():
    def __init__(self, model, inputsize, dataloader, lr=0.0001):
        self.model = model

        self.dataloader = dataloader
        self.trainloader = dataloader.trainloader
        self.valloader = dataloader.valloader

        self.bin_means = dataloader.bin_means

        #self.loss = self.calcLoss()
        #self.loss = torch.nn.MSELoss()
        #self.loss = torch.nn.L1Loss()
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)
        #self.optimizer = torch.optim.SGD(self.model.parameters(), lr=0.001)

        self.loss = self.calcLoss4
        self.accuracy = self.calcAccuracy2

        self.BCE = nn.BCELoss()
        #self.MSE = nn.MSELoss()
        #self.loss = nn.MSELoss()
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(self.device)


        self.model = self.model.to(self.device)
        print(model)
        self.total_params = sum(p.numel() for p in model.parameters())
        print("Total parameters: ", self.total_params)





    def train_one_epoch(self):
        #self.model.train()

        epoch_loss = 0
        for i, data in enumerate(self.trainloader):
            inputs, labels = data
            labels = labels[:, 3:]
            inputs = inputs.to(self.device)
            labels = labels.to(self.device)

            self.optimizer.zero_grad()

            preds, bases = self.model(inputs)

            loss = self.loss(preds, bases, labels)
            loss.backward()
            self.optimizer.step()

            epoch_loss += loss.item()

        avg_loss = (epoch_loss-loss.item()) / (len(self.trainloader)-1) # Remove the last batch, as it might not be full size
        return avg_loss


    def validate(self):
        self.model.train(False)
        #self.model.eval()

        loss_total = 0
        acc_total = 0
        for i, data in enumerate(self.valloader):
            inputs, labels = data
            labels = labels[:,3:]
            inputs = inputs.to(self.device)
            labels = labels.to(self.device)

            pred, base = self.model(inputs)

            #loss = self.calcLoss(pred, base, labels)
            loss = self.loss(pred, base, labels)
            loss_total += loss.item()
            acc = self.accuracy(pred, base, labels)

            acc_total += acc.item()

        loss_avg = (loss_total-loss.item()) / (len(self.valloader)-1)  # Remove the last batch, as it might not be full size
        acc_avg = (acc_total - acc.item()) / (len(self.valloader) - 1)  # Remove the last batch, as it might not be full size
        self.model.train(True)

        return loss_avg, acc_avg

    def calcEuclideanError(self, predictions, labels):
        errors = predictions.sub(labels)
        # error = predictions.sub(true_dF)
        sq_errors = torch.square(errors)
        sum_sq_errors = torch.sum(sq_errors, 1)
        sum_errors = torch.sqrt(sum_sq_errors)


        return sum_errors

    def calcLossScalars(self, base, labels):
        # Two ways of putting error in context of force
        true_dF = labels.sub(base)
        scalars = torch.sqrt(torch.sum(torch.square(true_dF), dim=1))  # Dunno which scalar to use, try both i guess
        #scalar = torch.div(scalar, len(labels))
        # scalar = torch.sqrt(torch.sum(torch.square(base)))
        scalars = torch.max(scalars, torch.tensor(0.00001))  # Dont let scalar become too small, avoid exploding loss
        return scalars

    def calcLoss1(self, predictions, base, labels):
        euclidean_errors = self.calcEuclideanError(predictions, labels)
        scalars = self.calcLossScalars(base, labels)
        scaled_errors = euclidean_errors.div(scalars)

        mean_err = torch.mean(scaled_errors)
        #mean_err = torch.median(scaled_errors)



        return mean_err

    def calcLoss2(self, predictions, base, labels):
        euclidean_errors = self.calcEuclideanError(predictions, labels)
        mean_err = torch.mean(euclidean_errors)

        return mean_err

    def calcLoss3(self, predictions, base, labels):
        return self.MSE(predictions, labels)

    def calcWeightedMean(self, x, weights):
        if torch.sum(weights == 0):
            print("Weights equal 0")
        return torch.mean(x*weights)/torch.sum(weights)

    def calcWeightedVariance(self, x, mean, weights):
        return torch.sum( torch.square(x-mean) * weights ) / torch.sum(weights)

    def calcVariances(self, pred, labels):
        labels = labels[:, 0]
        vars = []

        for i in range(pred.shape[1]):
            #mean = self.calcWeightedMean(labels, pred[:,i])
            mean = self.bin_means[i]
            var = self.calcWeightedVariance(labels, mean, pred[:,i])
            vars.append(var)

        vars = torch.stack(vars)
        vars = torch.nan_to_num_(vars)

        return vars


    def calcLoss4(self, pred, base, labels):
        return self.BCE(pred, labels)

        #variances = self.calcVariances(pred, labels)
        #return torch.sum(variances)

        #return torch.sum(torch.square(variances))


    def calcAccuracy(self, predictions, base, labels):
        euclidean_errors = self.calcEuclideanError(predictions, labels)
        scalars = self.calcLossScalars(base, labels)
        scaled_errors = euclidean_errors.div(scalars)


        ones = torch.ones(scaled_errors.shape, dtype=torch.float32).to(self.device)
        zeroes = torch.zeros(scaled_errors.shape, dtype=torch.float32).to(self.device)
        unbound_acc = ones.sub(scaled_errors)
        stack = torch.stack((unbound_acc, zeroes))
        bounded_acc = torch.max(stack, dim=0)[0]

        mean_acc = torch.mean(bounded_acc)
        return mean_acc

    def calcAccuracy2(self, pred, base, labels):
        return torch.tensor(0)

        #variances = self.calcVariances(pred, labels)
        #return torch.mean(variances)


    def saveModel(self, working_folder):
        path = working_folder + "model.pt"
        torch.save(self.model, path)


    def getBinnedForces(self, dataloader):
        labels_all = []
        predictions_all = []

        for i, data in enumerate(dataloader):
            inputs, labels = data
            labels = labels[:, 0]
            inputs = inputs.to(self.device)
            labels = labels.to(self.device)

            pred, base = self.model(inputs)

            labels_all.append(labels)
            predictions_all.append(pred)

        # labels_all = torch.cat(labels_all)[:,0]
        labels_all = torch.cat(labels_all)
        predictions_all = torch.cat(predictions_all, dim=0)

        predictions = torch.argmax(predictions_all, dim=1)

        forces = []
        for i in range(predictions_all.shape[1]):
            forces.append(labels_all[torch.where(predictions == i)].cpu().numpy())

        return forces

    def visualizePredictions(self):

        force_bins_training = self.getBinnedForces(self.trainloader)
        force_bins_validation = self.getBinnedForces(self.valloader)


        plt.figure()

        for i in range(len(force_bins_training)):
            plt.subplot(211)
            plt.hist(force_bins_training[i], bins=8, histtype='step')
            plt.yscale('log')
            plt.title("Training data")
            plt.subplot(212)
            plt.hist(force_bins_validation[i], bins=8, histtype='step')
            plt.yscale('log')
            plt.title("Validation data")


        plt.show()








"""
           #if (mean_err == 'inf'):
        #print(predictions)
        if (mean_err > 100000):
            print(mean_err)
            print(predictions)
            print(labels)
            print(euclidean_errors)
            exit()

        if (torch.sum(torch.isnan(mean_err))):
            print("Nan numbers!")
            print(mean_err)
            print(predictions)
            print(labels)
            exit()
"""