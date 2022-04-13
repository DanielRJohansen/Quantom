import numpy as np
import matplotlib.pyplot as plt
from time import time


class Trainer():
    def __init__(self, net, working_folder, n_neighbors):
        self.net = net
        self.working_folder = working_folder
        self.n_neighbors = n_neighbors

        self.best_acc = 0

    def train(self, n_epochs):
        val_loss_hist = []
        val_acc_hist = []
        train_loss_hist = []
        train_acc_hist = []


        for e in range(n_epochs):
            t0 = time()
            if (e == 0):                # Dont train the first epoch, so we can see a baseline!
                train_loss_hist.append(1)
            else:
                train_loss_hist.append(self.net.train_one_epoch())
            vloss, vacc = self.net.validate()
            val_loss_hist.append(vloss)
            val_acc_hist.append(vacc)


            if (vacc > self.best_acc):
                self.net.saveModel(self.working_folder)

            print("Epoch ", e, " train-loss: ", train_loss_hist[-1], " val-loss: ", val_loss_hist[-1], " mean accuracy: ", val_acc_hist[-1], " time: ", time()-t0)
        
        
        self.makePlot(train_loss_hist, train_acc_hist, val_loss_hist, val_acc_hist, n_epochs)
        
        

    def makePlot(self, train_loss, train_acc, val_loss, val_acc, n_epochs):
        train_loss[0] = train_loss[1]   # Not necessary, just looks better since we dont have data for epoch 0
        t = np.arange(0,n_epochs,1, dtype=np.int16)

        plt.plot(t, train_loss, "b--", t, val_loss, "r--", t, val_acc, "g-")
        plt.xlabel("Epoch")
        plt.ylabel("Accuracy/Loss")
        plt.legend(["train loss", "validation loss", "validation accuracy"])
        plt.title("Training - " + str(self.n_neighbors) + " neighbors - "
                  + str(self.net.total_params) + " parameters")
        plt.grid()
        plt.ylim(0,3)

        plt.savefig(self.working_folder + "\\train_plot_"+str(self.n_neighbors)+"N.png")
        plt.show()
        



    #def loadModel(self, path):


