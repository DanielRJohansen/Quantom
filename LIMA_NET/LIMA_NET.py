import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import numpy as np

from LIMA_TOOLS import *



class LIMANET():
    def __init__(self, model, inputsize, dataloader, bins=None, lr=0.0001):
        self.model = model

        self.dataloader = dataloader
        self.trainloader = dataloader.trainloader
        self.valloader = dataloader.valloader

        self.bins = bins


        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)
        #self.optimizer = torch.optim.SGD(self.model.parameters(), lr=lr)

        self.loss = self.calcLoss7
        self.accuracy = self.calcAccuracy

        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(self.device)

        #self.BCE = nn.BCELoss(weight=dataloader.bin_weights.to(self.device))
        #self.CELoss = nn.CrossEntropyLoss()
        self.CELoss = F.binary_cross_entropy_with_logits
        self.dir_vectors = getDirVectors().to(self.device)                                      # Vectors in a set of directions

        self.model = self.model.to(self.device)
        print(model)
        self.total_params = sum(p.numel() for p in model.parameters())
        print("Total parameters: ", self.total_params)





    def train_one_epoch(self):
        #self.model.train()

        epoch_loss = 0
        for i, data in enumerate(self.trainloader):
            inputs, labels = data

            #labels = labels[:, 3:]                                  # First 3 values are force xyz, next n values are onehot encodings
            #labels = labels[:,0]                                    # only forward x force

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
        acc_total = torch.zeros(2, dtype=float)
        #acc_total = 0
        for i, data in enumerate(self.valloader):
            inputs, labels = data

            inputs = inputs.to(self.device)
            labels = labels.to(self.device)

            pred, base = self.model(inputs[:,0,:,:])
            pred_baddata, _ = self.model(inputs[:,1,:,:])

            loss = self.loss(pred, base, labels)
            loss_total += loss.item()
            acc = self.accuracy(pred, pred_baddata, labels)

            acc_total += acc

        loss_avg = (loss_total-loss.item()) / (len(self.valloader)-1)  # Remove the last batch, as it might not be full size
        acc_avg = (acc_total - acc) / (len(self.valloader) - 1)  # Remove the last batch, as it might not be full size
        self.model.train(True)

        return loss_avg, acc_avg

    def calcEuclideanError(self, predictions, labels):
        errors = predictions.sub(labels)
        # error = predictions.sub(true_dF)
        sq_errors = torch.square(errors)
        sum_sq_errors = torch.sum(sq_errors, dim=1)
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


    def calcLoss6(self, pred, base, labels):
        dir_vector = self.generatePredVector(pred)
        label_vectors = normalizedVectors(labels[:, :3])

        dotp1 = dotP(dir_vector, label_vectors, dim=1)

        err = 1 - dotp1

        return torch.mean(err)

        dif = labels-pred
        confidence_undershot = torch.clamp(dif, min=0)
        confidence_overshot = torch.clamp(dif, max=0)

        err = torch.abs(confidence_overshot) + confidence_undershot#*confidence_undershot
        return torch.mean(err)
        #return torch.mean(dif)
        #dif[torch.where(dif, 0)] = 0
        print(dif.shape)
        exit()

        return self.CELoss(pred, labels)


    def calcLoss7(self, pred, base, labels):
        #print(pred.shape)
        #print(labels.shape)

        gauss_activations = self.calcGaussValue(pred[:,:,0], pred[:,:,1], labels[:,3:])
        #print("Gauss act: ", gauss_activations)
        #print(torch.min(gauss_activations))
        #log_activation = torch.log(gauss_activations)
        cost = 1 - gauss_activations/100


        mean_cost = torch.mean(cost)
        #print("Loss: ", mean_cost)
        #exit()
        return mean_cost



    def generatePredVector(self, pred):
        batch_size = pred.shape[0]
        summed_dir = torch.zeros((batch_size, 3)).to(self.device)

        for i in range(6):
            dir = self.dir_vectors[i]
            dir = dir.repeat(batch_size, 1)

            pred_confidence = pred[:,i].repeat(3,1).transpose(0,1)
            vector_contributions = pred_confidence * dir
            summed_dir += vector_contributions

        dir_vector = normalizedVectors(summed_dir).nan_to_num() # vectors of len 0 becoes (0,0,0)

        return dir_vector


    def calcAccuracy(self, pred, pred_baddata, labels):
        gauss_activations = self.calcGaussValue(pred[:,:,0], pred[:,:,1], labels[:,3:]) / 100

        gauss_activations_baddata = self.calcGaussValue(pred_baddata[:,:,0], pred_baddata[:,:,1], labels[:,3:]) / 100
        #acc = (1 / pred[:, :, 1]) * gauss_activations

        #gauss_activations = torch.clamp(gauss_activations, max=1)

        return torch.tensor((torch.mean(gauss_activations), torch.mean(gauss_activations_baddata)))







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

    def calcGaussValue(self, mean, std, x):

        #print(torch.min(std))
        # Variant of gauss with area=1
        exponent = -0.5 * ((x-mean)/std)**2
        activation = (1 / std) * torch.exp(exponent)
        #print(activation)
        return activation




    def noPred(self, dataloader):
        labels_all = []
        #predictions_all = []

        x_prev = []

        for i, data in enumerate(dataloader):
            inputs, labels = data

            labels_all.append(labels[:,0])
            x_prev.append(inputs[:,:,0])

        labels_all = torch.cat(labels_all, dim=0)
        x_prev = torch.cat(x_prev, dim=0)


        #print("xprev: ", x_prev.shape)
        onehots = onehotEncoder(x_prev, self.bins.bin_centers)



        forces = []
        for i in range(onehots.shape[1]):
            #forces.append(labels_all[torch.where(onehots == i)].cpu().numpy())
            forces.append(labels_all[torch.where(onehots[:,i])].cpu().numpy())

        return forces


        # labels_all = torch.cat(labels_all)[:,0]
        #labels_all = torch.cat(labels_all)
        #predictions_all = torch.cat(predictions_all, dim=0)

        #predictions = torch.argmax(predictions_all, dim=1)




    def doHists(self, force_bins_training, force_bins_validation):
        plt.figure()

        for i in range(len(force_bins_training)):
            plt.subplot(211)
            plt.hist(force_bins_training[i], bins=10, histtype='step')
            plt.yscale('log')
            plt.title("Training data")

            plt.subplot(212)
            plt.hist(force_bins_validation[i], bins=10, histtype='step')
            plt.yscale('log')
            plt.title("Validation data")



        plt.show()




    def visualizePredictions(self):
        self.model.train(False)
        # self.model.eval()

        loss_total = 0
        acc_total = torch.zeros(2, dtype=float)
        # acc_total = 0

        predictions = []
        labels = []

        for i, data in enumerate(self.valloader):
            inputs, label = data

            inputs = inputs.to(self.device)
            label = label.to(self.device)

            pred, base = self.model(inputs[:, 0, :, :])

            predictions.append(pred)
            labels.append(label[:,3:])



        predictions = torch.cat(predictions, dim=0)
        predictions = predictions.cpu().detach().numpy()

        labels = torch.cat(labels, dim=0)
        labels = labels.cpu().numpy()

        for i in range(6):
            plt.hist(labels[:, i], bins=20, alpha=0.3, label=str(i))
        plt.xlabel("True vector-activation mean")
        plt.ylabel("Occurances")
        plt.show()

        for i in range(6):
            plt.hist(predictions[:, i, 0], bins=20, alpha=0.3, label=str(i))
        plt.xlabel("Predicted vector-activation mean")
        plt.ylabel("Occurances")
        plt.show()

        for i in range(6):
            plt.hist(predictions[:, i, 1], bins=20, alpha=0.3, label=str(i))
        plt.ylabel("Occurances")
        plt.xlabel("Predicted vector-activation standard deviation")
        plt.show()

    def visualizePredictions1(self):

        force_bins_training = self.getBinnedForces(self.trainloader)
        force_bins_validation = self.getBinnedForces(self.valloader)
        forcebins_nopred = self.noPred(self.trainloader)



        self.doHists(force_bins_training, force_bins_validation)
        #self.doHists(force_bins_training, forcebins_nopred)

        plt.figure()
        plt.subplot(311)
        plt.boxplot(force_bins_training, vert=False)
        #print("mean ", np.mean(force_bins_training[i]))
        #print("std", np.std(force_bins_training[i]))
        #print("std", np.std(force_bins_validation[i]))
        plt.title("Training data")


        plt.subplot(312)
        plt.boxplot(force_bins_validation, vert=False)
        plt.title("Validation data")


        plt.subplot(313)
        plt.boxplot(forcebins_nopred, vert=False)
        plt.title("No prediction data")



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
            
            
            
            
            
            
            
                    #dir_vector = self.generatePredVector(pred)
        #dir_vector_baddata = self.generatePredVector(pred_baddata)
        #label_vectors = normalizedVectors(labels[:,:3])

        #euc_err = 1 - self.calcEuclideanError(dir_vector, label_vectors) * 0.5
        #euc_err_baddata = 1 - self.calcEuclideanError(dir_vector_baddata, label_vectors) * 0.5
        # return torch.tensor((torch.mean(euc_err), torch.mean(euc_err_baddata)))
        #print(dir_vector.shape)
        #exit()
        #dotp1 = dotP(dir_vector, label_vectors, dim=1)
        #dotp2 = dotP(dir_vector_baddata, label_vectors, dim=1)

        #print(torch.mean(dir_vector, dim=0))

        #return torch.tensor((torch.mean(dotp1), torch.mean(dotp2)))
"""