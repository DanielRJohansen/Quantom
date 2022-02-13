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
        raw_data = np.genfromtxt(data_filepath, delimiter=';', dtype=np.float32)

        labels = torch.from_numpy(raw_data[:, 0:3])
        data = torch.from_numpy(raw_data[:, 3:final_index])  # -1 removes the na that numpy adds for some odd reason..
        datapoints_total = data.shape[0]

        self.inputsize = data.shape[1]

        print("First label", labels[0,:])
        print("N input_floats", self.inputsize)
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


