from Dataloader import WaterforceDataloader
from LIMA_DNN import LIMADNN
from LIMA_NET import LIMANET
from torch.utils.data import DataLoader











if __name__ == '__main__':
    print("Hello world!")
    dataloader = WaterforceDataloader(batch_size=64)
    #dataloader = DataLoader(dataset=dataset, batch_size=4, shuffle=True, num_workers=0)

    model = LIMADNN(inputsize = dataloader.inputsize, dataloader=dataloader)
    print(model)
    #model.train()
    net = LIMANET(model, dataloader.inputsize, dataloader)
    net.train()