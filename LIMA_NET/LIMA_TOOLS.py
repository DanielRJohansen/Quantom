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