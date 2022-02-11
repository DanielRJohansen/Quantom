from LIMA_NET.Dataloader import WaterforceDataloader
from LIMA_DNN import LIMADNN
from torch.utils.data import DataLoader











if __name__ == '__main__':
    print("Hello world!")
    dataloader = WaterforceDataloader(batch_size=16)
    #dataloader = DataLoader(dataset=dataset, batch_size=4, shuffle=True, num_workers=0)

    net = LIMADNN(inputsize = dataloader.inputsize, dataloader=dataloader)
    print(net)
    net.train()

