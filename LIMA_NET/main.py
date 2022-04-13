from Dataloader import WaterforceDataloader
from LIMA_DNN import *
from LIMA_NET import LIMANET
from LIMA_TRAINER import  Trainer
from torch.utils.data import DataLoader
import torch










if __name__ == '__main__':
    n_neighbors = 8

    workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_200000\\"
    data_filepath = workdir + "traindata.bin"

    neighbors_per_row = 128

    net_type = 2

    dataloader = WaterforceDataloader(data_filepath, neighbors_per_row,batch_size=32, nearest_n_atoms=n_neighbors, prepper_id=net_type)

    if (net_type == 1):
        model = LIMADNN1(n_neighbors=n_neighbors)
    if (net_type == 2):
        model = LIMADNN2(n_neighbors=n_neighbors)


    net = LIMANET(model, dataloader.inputsize, dataloader)

    trainer = Trainer(net, workdir, n_neighbors=n_neighbors)

    trainer.train(200)