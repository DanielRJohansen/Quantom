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




def prepData(input_data, nearest_n_atoms):
    self_data = input_data[:,3:9]                                               # 3 force-prev, 3 deltaforce-prev
    self_data = self_data.reshape((self_data.shape[0], self_data.shape[1]//3, 3))

    neighbor_data = input_data[:, 9:]  # Remove query atom on all steps/rows  9 values because label has been removed!
    neighbor_data = neighbor_data.reshape((neighbor_data.shape[0], neighbor_data.shape[1] // 12, 12))       # Place neighors along new dimension
    neighbor_data = neighbor_data[:, :, :3]  # Remove all neighbor data except pos

    data = np.concatenate((self_data, neighbor_data), axis=1)           # First 2 units of dim 1 is self_particle : forceprev, dforceprev

    #print(data[:,:2,:])

    return torch.from_numpy(data)





class WaterforceDataloader():
    def __init__(self, data_filepath, neighbors_per_row, nearest_n_atoms=128, batch_size=64, norm_data=True, aug_data=False):


        values_per_row = 3 * 4 * (1 + neighbors_per_row)        # 2 for self pos + self_pos_prev

        raw_data = np.fromfile(data_filepath, dtype=np.float32)
        raw_data = raw_data.reshape((int(len(raw_data)/(values_per_row)), values_per_row))          # Reshape into steps
        raw_data = raw_data[:,0:3 * 4 * (1 + nearest_n_atoms)]                                              # Remove distant atoms


        # t0 forces x, y, z is always the label, which is why 0:3 below
        labels = torch.from_numpy(raw_data[:,0:3])
        input_data = raw_data[:,3:]

        data = prepData(input_data, nearest_n_atoms)



        datapoints_total = data.shape[0]
        self.inputsize = data.shape[1]*data.shape[2]
        self.printMetaInfo(data, labels, datapoints_total)
        self.n_train = int(datapoints_total * 0.8)


        directionvectors_labels = self.makeDiscreteDirectionsLabels(labels)

        labels = torch.cat((labels, directionvectors_labels), dim=1)

        #if norm_data:          # Doesn't work right now...
        #    data = self.normalizeData(data, self.n_train)                   # Normalize AFTER making labels



        #### Now split the dataset
        train_data = data[:self.n_train, :, :]
        train_labels = labels[:self.n_train, :]
        train_data_aug, train_labels_aug = self.generateFakeData(train_data, train_labels)
        train_data = torch.cat((train_data, train_data_aug), dim=0)
        train_labels = torch.cat((train_labels, train_labels_aug), dim=0)
        trainset = LIMADataset(train_data, train_labels)
        #trainset = LIMADataset(data[:self.n_train, :], labels[:self.n_train, :])



        val_data = data[self.n_train:, :, :]
        val_labels = labels[self.n_train:, :]
        val_data_obscured = self.obscureNeighbordata(val_data)          # Copy of the val_data but all neighbor position data is far away and thus useless
        val_data = torch.cat((  val_data[:, None, :, :], val_data_obscured[:, None, :, :]), dim=1)

        valset = LIMADataset(val_data, val_labels)
        #valset_obscured = LIMADataset(val_data_obscured, val_labels)

        self.batch_size = batch_size
        self.trainloader = DataLoader(dataset=trainset, batch_size=batch_size, shuffle=True, num_workers=2)
        self.valloader = DataLoader(dataset=valset, batch_size=batch_size, shuffle=False, num_workers=2)
        #self.valloader_obscured = DataLoader(dataset=valset_obscured, batch_size=batch_size, shuffle=False, num_workers=2)
        print("\n\n")




        print("Inputdata shape ", data.shape)
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




    def makeDiscreteDirectionsLabels(self, labels):         #(elements, 3)
        dir_vectors = getDirVectors()
        n_vectors = dir_vectors.shape[0]
        n_datapoints = labels.shape[0]



        #print(labels.shape)
        #print(dir_vectors.shape)

        dirlabels = []

        force_normalized = normalizedVectors(labels)



        for i in range(n_vectors):
            dir_vector = dir_vectors[i,:].repeat(n_datapoints, 1)#.transpose(0,1)   # Select dir vector, repeat for each datapoint

            #force_dir_likeness = torch.sum(dir_vector * force_normalized, dim=1)        # Dotproduct from Stackoverflow lol
            force_dir_likeness = dotP(dir_vector, force_normalized, dim=1)  # Dotproduct from Stackoverflow lol
            force_dir_likeness = torch.maximum(force_dir_likeness, torch.tensor(0))     # This allows the network to not be punished for pushing from both sides!
            dirlabels.append(force_dir_likeness)

        dirlabels = torch.stack(dirlabels).transpose(0,1)
        print("Direction labels shape: ", dirlabels.shape)
        dirlabels = F.softmax(dirlabels, dim=1)                                         # We need to softmax this, because the NN also softmaxes its output. They must be able to match
        #exit()
        return dirlabels




    def generateVectorDistributionFromSample(self, sample):
        #print("Sample: ", sample.shape)
        n_datapoints = sample.shape[0]
        data_force_scalars = vectorLengths(sample[:, 0, :])
        force_mean = torch.mean(data_force_scalars).repeat(n_datapoints)
        force_std = torch.std(data_force_scalars).repeat(n_datapoints)
        fake_force_scalars = torch.normal(force_mean, force_std)
        fake_force_scalars = fake_force_scalars.repeat(3, 1).transpose(0, 1)
        fake_force_dirs = normalizedVectors(sampleNormalizedData((n_datapoints, 3), 0, 1))
        fake_forces = fake_force_dirs * fake_force_scalars
        print("Mean, Std: ", force_mean[0], force_std[0])
        return fake_forces


    def getFakePositions(self, dims):
        fake_positions = sampleNormalizedData(dims, mean=0, std=1)
        fake_positions = normalizeVectors2(fake_positions) * 4  # So they are all a distance of 4 nm away from p0
        fake_positions = fake_positions + sampleNormalizedData(fake_positions.shape, mean=0, std=0.5)
        return fake_positions


    def generateFakeData(self, data, labels, fakedata_ratio=0.2, preserve_forces=False):
        n_datapoints = data.shape[0]
        n_neighbors = data.shape[1]-2   # First two is p0 data
        n_select = int(n_datapoints*fakedata_ratio)


        fake_forces_data = self.generateVectorDistributionFromSample(data)

        fake_positions = self.getFakePositions((n_datapoints, n_neighbors, 3))

        new_data = torch.cat((fake_forces_data[:, None, :], fake_forces_data[:, None, :], fake_positions), dim=1) # Put forces twice since second element is unused anyways
        new_data = new_data[:n_select, :, :]




        fake_forces_labels = self.generateVectorDistributionFromSample(labels[:, None, :3])
        fake_forces_labels_dv = self.makeDiscreteDirectionsLabels(fake_forces_labels)
        new_labels = torch.cat((fake_forces_labels, fake_forces_labels_dv), dim=1)
        new_labels = new_labels[:n_select, :]

        return new_data, new_labels

    def obscureNeighbordata(self, data):
        n_datapoints = data.shape[0]
        n_neighbors = data.shape[1] - 2  # First two is p0 data

        fake_positions = self.getFakePositions((n_datapoints, n_neighbors, 3))

        print(fake_positions.shape)
        untouched_data = data[:,:2,:]
        print(untouched_data.shape)

        new_data = torch.cat((untouched_data, fake_positions), dim=1)
        return new_data




    def assignWeightsToBins(self, onehot_labels):
        #print(onehot_labels.shape)

        bin_counts = torch.sum(onehot_labels[:self.n_train, :], dim=0)

        bin_weights = torch.tensor(1e+6, dtype=float) / bin_counts
        bin_weights = torch.nan_to_num_(bin_weights, 1, 1, 1)  # turns nan inf -inf to 1
        #bin_weights = torch.nn.Softmax(bin_weights)
        print("Bin weights: ", bin_weights)
        return bin_weights



    def normalizeData(self, data, n_train):                         # Needs to be fixed, after adding force back in data
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