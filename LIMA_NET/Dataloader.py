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



def prepData1(input_data, nearest_n_atoms):
    #prev_data = raw_data[:, 3:]
    neighbor_data = input_data[:, 3:]
    #       labels = torch.from_numpy(raw_data[:, 0:3])
    #      data = torch.from_numpy(raw_data[:, 3:final_index])  # -1 removes the na that numpy adds for some odd reason..

    # All this does is rearrange data so all postions of neighbors comes first, then all forces of neighbors
    neighbor_data = neighbor_data.reshape((neighbor_data.shape[0], neighbor_data.shape[1] // 3, 3))
    pos_data = neighbor_data[:, ::2, :].reshape((neighbor_data.shape[0], nearest_n_atoms * 3))
    force_data = neighbor_data[:, 1::2, :].reshape((neighbor_data.shape[0], nearest_n_atoms * 3))

    data = np.concatenate((input_data, pos_data, force_data), axis=1)
    return torch.from_numpy(data)

def prepData2(input_data, nearest_n_atoms):
    ones = np.ones((input_data.shape[0], 3), dtype=np.float32)      # Insert ones where the label was
    print(input_data.shape)
    data = np.concatenate((ones, input_data), axis=1)
    print(data.shape)
    data = data.reshape((data.shape[0], data.shape[1] // 12, 12))     #
    print("First row ", data[0,0,:])
    print(data.shape)
    exit()
    return torch.from_numpy(data)

def prepData3(input_data, nearest_n_atoms):

    return torch.from_numpy(input_data[:, 3:])


def prepData4(input_data, nearest_n_atoms):
    data = input_data[:,9:]                                             # Remove query atom on all steps/rows
    data = data.reshape((data.shape[0], data.shape[1] // 12, 12))       # Place neighors along new dimension
    data = data[:,:,:3]                                                 # Remove all neighbor data except pos
    print("Data shape: ",  data.shape)
    return torch.from_numpy(data)



class WaterforceDataloader():
    def __init__(self, data_filepath, neighbors_per_row, nearest_n_atoms=128, batch_size=64, prepper_id=0, bins=16, n_queries=1):
        preppers = [prepData1, prepData2, prepData3, prepData4]
        data_prepper = preppers[prepper_id-1]

        values_per_row = 3 * 4 * (1 + neighbors_per_row)        # 2 for self pos + self_pos_prev

        raw_data = np.fromfile(data_filepath, dtype=np.float32)
        #print("Max value in dataset ", max(raw_data))
        #print("Min value in dataset ", min(raw_data))


        raw_data = raw_data.reshape((int(len(raw_data)/(values_per_row)), values_per_row))          # Reshape into steps
        raw_data = raw_data[:,0:3 * 4 * (1 + nearest_n_atoms)]                                              # Remove distant atoms

        labels = torch.from_numpy(raw_data[:,0:3])
        input_data = raw_data[:,3:]

        data = data_prepper(input_data, nearest_n_atoms)
        #data = torch.from_numpy(data)



        datapoints_total = data.shape[0]

        self.inputsize = data.shape[1]

        print("Input data shape ", data.shape)
        print("First label ", labels[0,:])
        print("First base-force ", data[0,0:3])
        print("N input values", self.inputsize)
        print("Final input value check", data[0, -1])
        print("N datapoints total: ", datapoints_total)










        #### Now split the dataset
        self.n_train = int(datapoints_total * 0.85)






        # Now make bins of the data
        self.datapoint_max = torch.max(torch.abs(labels[0:self.n_train, 0]))     # temp, now it only checks x forces!
        self.datapoint_min = self.datapoint_max*-1
        self.bin_means = self.calcBinMeans(bins)
        onehot_labels = self.makeBinaryLabels(labels, self.bin_means)


        #labels = [labels, onehot_labels]
        labels = torch.cat((labels, onehot_labels), dim=1)    # First 3 values are force xyz, next n values are onehot encodings
        #print(labels.shape)
        #exit()



        trainset = LIMADataset(data[0:self.n_train, :], labels[0:self.n_train, :])
        valset = LIMADataset(data[self.n_train:, :], labels[self.n_train:, :])
        self.batch_size = batch_size
        self.trainloader = DataLoader(dataset=trainset, batch_size=batch_size, shuffle=True, num_workers=2)
        self.valloader = DataLoader(dataset=valset, batch_size=batch_size, shuffle=False, num_workers=2)
        print("\n\n")



    def calcBinMeans(self, n_bins):
        min = self.datapoint_min
        max = self.datapoint_max

        bin_stride = (max-min)/(n_bins+1)
        bin_means = []
        for i in range(n_bins):
            bin_means.append(min + i * bin_stride + bin_stride * 0.5)
        print(bin_means)
        return bin_means

    def makeBinaryLabels(self, labels, bin_means):
        #errors = torch.zeros((labels.shape[0], len(bin_means)), dtype=torch.float32)
        errors = []

        for i in range(len(bin_means)):
            errors.append(torch.square(labels[:,0]-bin_means[i]))
        errors = torch.stack(errors, 1)

        max_vals = torch.min(errors, dim=1).values.view(-1,1)
        max_vals = torch.repeat_interleave(max_vals, len(bin_means), dim=1)

        labels_ = torch.where(errors==max_vals, True, False)
        print(labels_)
        print("a shape ", labels_.shape)
        print(torch.sum(labels_, dim=0))

        return labels_.float()
