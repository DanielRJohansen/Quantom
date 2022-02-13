from Dataloader import WaterforceDataloader
from LIMA_DNN import LIMADNN
from LIMA_NET import LIMANET
from torch.utils.data import DataLoader
import torch










if __name__ == '__main__':
    n_neighbors = 1
    data_filepath = "D:\\Quantom\\LIMANET\\sim_out\\atom0_lines29556.csv"


    dataloader = WaterforceDataloader(data_filepath, batch_size=32, nearest_n_atoms=n_neighbors)

    model = LIMADNN(n_neighbors=n_neighbors)
    print(model)
    pytorch_total_params = sum(p.numel() for p in model.parameters())
    print("Total parameters: ", pytorch_total_params)
    #model.train()
    net = LIMANET(model, dataloader.inputsize, dataloader)
    net.train(5000)