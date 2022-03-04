from Dataloader import WaterforceDataloader
from LIMA_DNN import LIMADNN
from LIMA_NET import LIMANET
from LIMA_TRAINER import  Trainer
from torch.utils.data import DataLoader
import torch










if __name__ == '__main__':
    n_neighbors = 0

    workdir = "C:\\PROJECTS\\Quantom\\LIMANET\\sim_out\\Steps_10000"
    data_filepath = workdir + "\\traindata.bin"

    dataloader = WaterforceDataloader(data_filepath, batch_size=32, nearest_n_atoms=n_neighbors)

    model = LIMADNN(n_neighbors=n_neighbors)

    #model.train()
    net = LIMANET(model, dataloader.inputsize, dataloader)

    trainer = Trainer(net, workdir, n_neighbors=n_neighbors)

    trainer.train(100)
    #net.train(5000)