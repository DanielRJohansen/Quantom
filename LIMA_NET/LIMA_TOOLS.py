import torch



def onehotEncoder(values, bin_means):
    # errors = torch.zeros((labels.shape[0], len(bin_means)), dtype=torch.float32)
    errors = []

    print("values ", values.shape)
    print("binmeans ", bin_means.shape)
    for i in range(len(bin_means)):
        errors.append(torch.square(values[:, 0] - bin_means[i]))
    errors = torch.stack(errors, 1)

    max_vals = torch.min(errors, dim=1).values.view(-1, 1)
    max_vals = torch.repeat_interleave(max_vals, len(bin_means), dim=1)

    labels_ = torch.where(errors == max_vals, True, False)
    return labels_



def getDirVectors():
    directions = torch.zeros((6, 3))   # 6 vectors with 3 dims
    for dim in range(3):
        for i in range(2):
            directions[i + dim*2, dim] = -1 + 2*i   # alternate sign for consecutive vectors: [-1,0,0], [1,0,0], [0,-1..
    return directions

def vectorLengths(array, dim=1):    # (n_vectors, n_dims(3))
    sq = torch.square(array)
    sums = torch.sum(sq, dim=dim)
    lengths = torch.sqrt(sums)  # We do NOT check if any forces are 0, this would prob be an error somewhere elsooo
    return lengths

def normalizedVectors(array):
    lengths = vectorLengths(array)# We do NOT check if any forces are 0, this would prob be an error somewhere elsooo
    #print(array[-20:,:])
    #print(lengths[-20:])
    #print(lengths.shape)
    #print(torch.count_nonzero(lengths))
    lengths = lengths.repeat(3, 1).transpose(0, 1)
    normalized_vectors = array/lengths
    #print(normalized_vectors)
    #exit()
    return normalized_vectors

def normalizeVectors2(array):
    #print(array.shape)
    lengths = vectorLengths(array, dim=2)# We do NOT check if any forces are 0, this would prob be an error somewhere elsooo


    lengths = lengths.repeat(3, 1, 1).transpose(0, 1).transpose(1,2)
    normalized_vectors = array/lengths

    return normalized_vectors


def sampleNormalizedData(dims, mean, std):
    mean = torch.ones(dims) * mean
    std = torch.ones(dims) * std
    sample = torch.normal(mean, std)
    return sample

def dotP(v1, v2, dim):                  # Dotproduct from Stackoverflow lol
    return torch.sum(v1 * v2, dim=dim)