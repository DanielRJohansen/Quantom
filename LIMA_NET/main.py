from Dataloader import WaterforceDataloader
from LIMA_DNN import LIMADNN
from LIMA_NET import LIMANET
from torch.utils.data import DataLoader
import torch










if __name__ == '__main__':
    dataloader = WaterforceDataloader(batch_size=32, nearest_n_atoms=16)

    model = LIMADNN(n_neighbors=4)
    print(model)
    pytorch_total_params = sum(p.numel() for p in model.parameters())
    print("Total parameters: ", pytorch_total_params)
    #model.train()
    net = LIMANET(model, dataloader.inputsize, dataloader)
    net.train(500)