from Dataloader import WaterforceDataloader
from LIMA_DNN import LIMADNN
from LIMA_NET import LIMANET
from LIMA_TRAINER import  Trainer
from torch.utils.data import DataLoader
import torch










if __name__ == '__main__':
    n_neighbors = 2

    workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_150000\\"
    data_filepath = workdir + "traindata.bin"

    neighbors_per_row = 128
    dataloader = WaterforceDataloader(data_filepath, neighbors_per_row,batch_size=32, nearest_n_atoms=n_neighbors)

    model = LIMADNN(n_neighbors=n_neighbors)


    net = LIMANET(model, dataloader.inputsize, dataloader)

    trainer = Trainer(net, workdir, n_neighbors=n_neighbors)

    trainer.train(50)