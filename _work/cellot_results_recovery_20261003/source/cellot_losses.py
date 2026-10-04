import torch

def compute_loss_g(f, g, source, transport=None):
    if transport is None:
        transport = g.transport(source)

    return f(transport) - torch.multiply(source, transport).sum(-1, keepdim=True)

def compute_loss_f(f, g, source, target, transport=None):
    if transport is None:
        transport = g.transport(source)

    return -f(transport) + f(target)
