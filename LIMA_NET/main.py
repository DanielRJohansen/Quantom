from Dataloader import WaterforceDataloader
from LIMA_DNN import LIMADNN
from LIMA_NET import LIMANET
from LIMA_TRAINER import  Trainer
from torch.utils.data import DataLoader
import torch










if __name__ == '__main__':
    n_neighbors = 4
    #data_filepath = "D:\\Quantom\\LIMANET\\sim_out\\atom0_lines4865.csv"
    #data_filepath = "D:\\Quantom\\LIMANET\\sim_out\\atom0_lines4933_shuffled.bin"
    #data_filepath = "D:\\Quantom\\LIMANET\\sim_out\\atom0_lines49741_shuffled.bin"
    data_filepath = "D:\\Quantom\\LIMANET\\sim_out\\atom0_lines49741.bin"

    dataloader = WaterforceDataloader(data_filepath, batch_size=32, nearest_n_atoms=n_neighbors)

    model = LIMADNN(n_neighbors=n_neighbors)

    #model.train()
    net = LIMANET(model, dataloader.inputsize, dataloader)

    trainer = Trainer(net, "D:\\Quantom\\LIMANET\\training")

    trainer.train(40)
    #net.train(5000)