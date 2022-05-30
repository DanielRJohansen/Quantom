import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import numpy as np
import math

from LIMA_TOOLS import *

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

    #data = input_data[:, :9].reshape((input_data.shape[0], 1, 9))

    data = data[:,:,:3]                                                 # Remove all neighbor data except pos
    #data = np.concatenate((data[:,:,:3], data[:,:,6:9]), axis=2)                                         # Remove all neighbor data except pos

    #data = data[:,:,3:4]

    return torch.from_numpy(data)





class WaterforceDataloader():
    def __init__(self, data_filepath, neighbors_per_row, bin_centers=None, nearest_n_atoms=128, batch_size=64, prepper_id=0, norm_data=True):
        preppers = [prepData1, prepData2, prepData3, prepData4]
        data_prepper = preppers[prepper_id-1]

        values_per_row = 3 * 4 * (1 + neighbors_per_row)        # 2 for self pos + self_pos_prev

        raw_data = np.fromfile(data_filepath, dtype=np.float32)
        raw_data = raw_data.reshape((int(len(raw_data)/(values_per_row)), values_per_row))          # Reshape into steps
        raw_data = raw_data[:,0:3 * 4 * (1 + nearest_n_atoms)]                                              # Remove distant atoms


        # t0 forces x, y, z is always the label, which is why 0:3 below
        labels = torch.from_numpy(raw_data[:,0:3])
        input_data = raw_data[:,3:]

        data = data_prepper(input_data, nearest_n_atoms)
        #data = torch.from_numpy(data)



        datapoints_total = data.shape[0]
        self.inputsize = data.shape[1]*data.shape[2]


        self.printMetaInfo(data, labels, datapoints_total)



        self.n_train = int(datapoints_total * 0.8)


        # Now make bins of the data
        self.datapoint_max = torch.max(torch.abs(labels[0:self.n_train, 0]))     # temp, now it only checks x forces!
        self.datapoint_min = self.datapoint_max*-1
        #onehot_labels = self.makeBinaryLabels(labels, bin_centers)
        directionvectors_labels = self.makeDiscreteDirectionsLabels(labels)


        #self.bin_weights = self.assignWeightsToBins(onehot_labels)                  # The weights depend on the (training) dataset

        #labels = [labels, onehot_labels]
        #labels = torch.cat((labels, onehot_labels), dim=1)    # First 3 values are force xyz, next n values are onehot encodings
        labels = torch.cat((labels, directionvectors_labels), dim=1)

        if norm_data:
            data = self.normalizeData(data, self.n_train)                   # Normalize AFTER making labels

        #### Now split the dataset
        trainset = LIMADataset(data[:self.n_train, :], labels[:self.n_train, :])
        valset = LIMADataset(data[self.n_train:, :], labels[self.n_train:, :])
        self.batch_size = batch_size
        self.trainloader = DataLoader(dataset=trainset, batch_size=batch_size, shuffle=True, num_workers=2)
        self.valloader = DataLoader(dataset=valset, batch_size=batch_size, shuffle=False, num_workers=2)
        print("\n\n")




        print("Inputdata shape ", data.shape)

        #ts1 = data[:,0,0]
        ts1 = raw_data[:, 6]
        t0 = labels[:,0]

        dif = torch.abs(t0 - ts1)

        index = torch.argmax(dif)

        print("Max label dif. Index: ", index, ts1[index], t0[index], " Dif ", dif[index])



    def printMetaInfo(self, data, labels, datapoints_total):
        print("Input data shape ", data.shape)
        print("First label ", labels[0, :])
        print("First base-force ", data[0, 0:3])
        print("N input values", self.inputsize)
        print("Final input value check", data[0, -1])
        print("N datapoints total: ", datapoints_total)

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
        labels_ = onehotEncoder(labels, bin_means)
        print("Bin means: ", bin_means)
        print("Bin counts: ", torch.sum(labels_, dim=0))



        return labels_.float()


    def makeDiscreteDirectionsLabels(self, labels):
        dir_vectors = getDirVectors()
        n_vectors = dir_vectors.shape[0]
        n_datapoints = labels.shape[0]


        #print(labels.shape)
        #print(dir_vectors.shape)

        dirlabels = []

        force_normalized = normalizedVectors(labels)


        for i in range(n_vectors):
            dir_vector = dir_vectors[i,:].repeat(n_datapoints, 1)#.transpose(0,1)   # Select dir vector, repeat for each datapoint

            force_dir_likeness = torch.sum(dir_vector * force_normalized, dim=1)        # Dotproduct from Stackoverflow lol
            force_dir_likeness = torch.maximum(force_dir_likeness, torch.tensor(0))     # This allows the network to not be punished for pushing from both sides!
            dirlabels.append(force_dir_likeness)

        dirlabels = torch.stack(dirlabels).transpose(0,1)
        print("Direction labels shape: ", dirlabels.shape)
        dirlabels = F.softmax(dirlabels, dim=1)                                         # We need to softmax this, because the NN also softmaxes its output. They must be able to match
        #exit()
        return dirlabels






    def assignWeightsToBins(self, onehot_labels):
        #print(onehot_labels.shape)

        bin_counts = torch.sum(onehot_labels[:self.n_train, :], dim=0)

        bin_weights = torch.tensor(1e+6, dtype=float) / bin_counts
        bin_weights = torch.nan_to_num_(bin_weights, 1, 1, 1)  # turns nan inf -inf to 1
        #bin_weights = torch.nn.Softmax(bin_weights)
        print("Bin weights: ", bin_weights)
        return bin_weights



    def normalizeData(self, data, n_train):
        # Only use means/stds of train as we dont know of the val ;)
        # Each column (pos x, py, pz, force x, fy and so are normalized separately

        means = torch.mean(data[:n_train, :, :], dim=(0,1))     # x, y, z
        vars = torch.std(data[:n_train, :, :], dim=(0,1))

        print(data.shape)
        print(data[0,0,:])
        data = (data - means) / vars
        print(data.shape)
        print(data[0, 0, :])

        #exit()
        return data