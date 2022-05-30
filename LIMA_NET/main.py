from Dataloader import WaterforceDataloader
from LIMA_DNN import *
from LIMA_NET import LIMANET
from LIMA_TRAINER import  Trainer
from torch.utils.data import DataLoader
import torch

import numpy as np
from scipy.stats import norm


class Bins:
    def __init__(self, bins, bincenters):
        self.bins = bins
        self.bin_centers = bincenters


def getBins(n=7, min=-4e+5, max=4e+5):

    if not (n % 2):
        print("Bins must be uneven")
        exit()

    n = (n+1) // 2                      # Make binstart-stops for np histograms
    l = np.geomspace(-512, min, n)
    r = np.geomspace(512, max, n)

    bins = np.concatenate((np.flip(l), r))
    bincenters = 0.5 * (bins[1:] + bins[:-1])


    return Bins(bins, bincenters)






if __name__ == '__main__':

    n_bins = 7                          # Must be uneven
    bins = getBins(n_bins)
    normalize_data = True               # Should we norm labels aswell?



    workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_350000\\"
    #data_filepath = workdir + "traindata.bin"
    data_filepath = workdir + "traindata_queries30.bin"


    neighbors_per_row = 32                  # in dataset
    n_neighbors = 32                        # To train with - obviously cannot be > n per row

    #dataloader = WaterforceDataloader(data_filepath, neighbors_per_row,batch_size=256, nearest_n_atoms=n_neighbors, prepper_id=2)
    dataloader = WaterforceDataloader(data_filepath, neighbors_per_row, bin_centers=bins.bin_centers, batch_size=4096, nearest_n_atoms=n_neighbors, prepper_id=4, norm_data=normalize_data)

#    model = LIMADNN1(n_neighbors=n_neighbors)
 #   model = LIMADNN2(n_neighbors=n_neighbors)
    #model = LIMADNN3(n_neighbors=n_neighbors)
    model = LIMADNN4(n_neighbors=n_neighbors, n_bins=n_bins, inputsize=dataloader.inputsize)

    net = LIMANET(model, dataloader.inputsize, dataloader, bins, lr=3e-4)

    trainer = Trainer(net, workdir, n_neighbors=n_neighbors)

    trainer.train(15)