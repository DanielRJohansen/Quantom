import numpy as np
import matplotlib.pyplot as plt

class Trainer():
    def __init__(self, net, working_folder):
        self.net = net
        self.working_folder = working_folder



    def train(self, n_epochs):
        val_loss_hist = []
        val_acc_hist = []
        train_loss_hist = []
        train_acc_hist = []


        for e in range(n_epochs):
            if (e == 0):                # Dont train the first epoch, so we can see a baseline!
                train_loss_hist.append(1)
            else:
                train_loss_hist.append(self.net.train_one_epoch())
            vloss, vacc = self.net.validate()
            val_loss_hist.append(vloss)
            val_acc_hist.append(vacc)

            print("Epoch ", e, " train-loss: ", train_loss_hist[-1], " val-loss: ", val_loss_hist[-1], " mean accuracy: ", val_acc_hist[-1])
        
        
        self.makePlot(train_loss_hist, train_acc_hist, val_loss_hist, val_acc_hist, n_epochs)
        
        

    def makePlot(self, train_loss, train_acc, val_loss, val_acc, n_epochs):
        train_loss[0] = train_loss[1]   # Not necessary, just looks better since we dont have data for epoch 0
        t = np.arange(0,n_epochs,1, dtype=np.int16)

        plt.plot(t, train_loss, "b--", t, val_loss, "r--", t, val_acc, "g-")
        plt.xlabel("Epoch")
        plt.ylabel("Accuracy/Loss")
        plt.legend(["train loss", "validation loss", "validation accuracy"])
        plt.title("Training - " + str(self.net.total_params) + " parameters")
        plt.grid()

        plt.savefig(self.working_folder + "\\training.png")
        plt.show()
        
        

