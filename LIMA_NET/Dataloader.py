import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import math



class LIMADataset(Dataset):
    def __init__(self, data, labels):
        self.data = data
        self.labels = labels
        self.n_samples = self.data.shape[0]
        print("N datapoints ", self.n_samples)

    def __getitem__(self, index):
        return self.data[index], self.labels[index]

    def __len__(self):
        return self.n_samples


class WaterforceDataloader():
    def __init__(self, data_filepath, nearest_n_atoms=69, batch_size=64):
        final_index = 6 + 3 * nearest_n_atoms * 2
        values_per_row = 70*3*2

        raw_data = np.fromfile(data_filepath, dtype=np.double)
        print(raw_data[0:6])
        print(raw_data.shape)
        raw_data = raw_data.reshape((int(len(raw_data)/(values_per_row)), values_per_row))
        print(raw_data.shape)
        print(raw_data[0, 0:6])
        #raw_data = raw_data.transpose()


        raw_data = raw_data.astype(np.float32)

        #raw_data = np.genfromtxt(data_filepath, delimiter=';', dtype=np.float32)
        #print(raw_data.shape)
        #exit()
        #print(raw_data[0,0:6])
        #exit(0)

        labels = torch.from_numpy(raw_data[:, 0:3])
        data = torch.from_numpy(raw_data[:, 3:final_index])  # -1 removes the na that numpy adds for some odd reason..
        datapoints_total = data.shape[0]

        self.inputsize = data.shape[1]

        print("First label ", labels[0,:])
        print("First base-force ", data[0,0:3])
        print("N input values", self.inputsize)
        print("Final input value check", data[0, -1])
        print("N datapoints total: ", datapoints_total)


        #### Now split the dataset
        self.n_train = int(datapoints_total * 0.7)
       # n_val = int(datapoints_total * 0.3)
        trainset = LIMADataset(data[0:self.n_train, :], labels[0:self.n_train, :])
        valset = LIMADataset(data[self.n_train:, :], labels[self.n_train:, :])
        self.batch_size = batch_size
        self.trainloader = DataLoader(dataset=trainset, batch_size=batch_size, shuffle=True, num_workers=2)
        self.valloader = DataLoader(dataset=valset, batch_size=batch_size, shuffle=False, num_workers=2)
        print("\n\n")


