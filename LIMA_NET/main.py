from Dataloader import WaterforceDataloader
from LIMA_DNN import *
from LIMA_NET import LIMANET
from LIMA_TRAINER import  Trainer
from torch.utils.data import DataLoader
import torch










if __name__ == '__main__':

    workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_350000\\"
    #data_filepath = workdir + "traindata.bin"
    data_filepath = workdir + "traindata_queries30.bin"


    neighbors_per_row = 32                  # in dataset

    n_neighbors = 16                        # To train with
    n_bins = 32

    #dataloader = WaterforceDataloader(data_filepath, neighbors_per_row,batch_size=256, nearest_n_atoms=n_neighbors, prepper_id=2)
    dataloader = WaterforceDataloader(data_filepath, neighbors_per_row, batch_size=512, nearest_n_atoms=n_neighbors, prepper_id=4, bins=n_bins)

#    model = LIMADNN1(n_neighbors=n_neighbors)
 #   model = LIMADNN2(n_neighbors=n_neighbors)
    #model = LIMADNN3(n_neighbors=n_neighbors)
    model = LIMADNN4(n_neighbors=n_neighbors, n_bins=n_bins)

    net = LIMANET(model, dataloader.inputsize, dataloader, lr=3e-4)

    trainer = Trainer(net, workdir, n_neighbors=n_neighbors)

    trainer.train(50)