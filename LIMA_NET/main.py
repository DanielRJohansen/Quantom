from Dataloader import WaterforceDataloader
from LIMA_DNN import *
from LIMA_NET import LIMANET
from LIMA_TRAINER import  Trainer
from torch.utils.data import DataLoader
import torch

import numpy as np
from scipy.stats import norm





if __name__ == '__main__':

    n_bins = 6                          # Must be uneven
    normalize_data = True               # Should we norm labels aswell?



    workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_350000\\"
    #data_filepath = workdir + "traindata.bin"
    #data_filepath = workdir + "traindata_queries2.bin" # something is wrong 50k forces are 0 in particle2
    data_filepath = workdir + "traindata.bin"


    neighbors_per_row = 32                  # in dataset
    n_neighbors = 32                        # To train with - obviously cannot be > n per row

    #dataloader = WaterforceDataloader(data_filepath, neighbors_per_row,batch_size=256, nearest_n_atoms=n_neighbors, prepper_id=2)
    dataloader = WaterforceDataloader(data_filepath, neighbors_per_row,  batch_size=128, nearest_n_atoms=n_neighbors, norm_data=normalize_data)

#    model = LIMADNN1(n_neighbors=n_neighbors)
 #   model = LIMADNN2(n_neighbors=n_neighbors)
    #model = LIMADNN3(n_neighbors=n_neighbors)
    model = LIMADNN4(n_neighbors=n_neighbors, n_bins=n_bins, inputsize=dataloader.inputsize)

    net = LIMANET(model, dataloader.inputsize, dataloader, lr=1e-4)

    trainer = Trainer(net, workdir, n_neighbors=n_neighbors)

    trainer.train(1)