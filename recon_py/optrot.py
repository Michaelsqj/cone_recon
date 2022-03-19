"""
The code uses pytorch for gradient descent optimization of in-plane rotation.
torchkbnufft is used for gridding and calculating weight
"""

import torch
import torch.nn as nn
import torchkbnufft as tkbn
import numpy as np
import os
if (torch.cuda.is_available()):
    device = torch.device("cuda")
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'
    print("Use GPU")
else:
    device = torch.device("cpu")
    print("Use CPU")
import scipy.io as io
from tqdm import tqdm

class OptTheta(nn.Module):
    def __init__(self, Grs, baseK, im_size):
        super().__init__()
        self.GrPRS = torch.tensor(Grs, dtype=torch.float32)     # (N,3)
        self.im_size = im_size
        self.baseK = baseK   # (M,3)
        self.len = baseK.shape[0]
        self.NLines = Grs.shape[0]
        self.pts = self.len * self.NLines
        self.baseKp = torch.tensor(np.reshape(baseK[:,0], (1, 1, -1), order='F'), dtype=torch.float32)
        self.baseKr = torch.tensor(np.reshape(baseK[:,1], (1, 1, -1), order='F'), dtype=torch.float32)
        self.baseKs = torch.tensor(np.reshape(baseK[:,2], (1, 1, -1), order='F'), dtype=torch.float32) # (1,1,M)
        self.theta = torch.rand(self.NLines, 1, requires_grad=True) * 2 * np.pi
    def forward(self):
        self.calc_slice()
        self.rot_baseK()
        dcomp = tkbn.calc_density_compensation_function(ktraj=self.ktraj.to(device), im_size=self.im_size)
        w = torch.squeeze(dcomp)
        loss = torch.sqrt(self.pts * torch.sum(w*w)) / torch.sum(w)
        return loss

    def calc_slice(self):
        """
        Args:
            GrPRS: (N,3)
            theta: (N,1)
        Returns:
            GpRad: (N,3)
            GrRad: (N,3)
            GsRad: (N,3)
        """
        Gr0 = torch.tensor([[-1., 0., 0.]])     # (1,3)
        Gs0 = torch.tensor([[0., 0., 1.]])
        iR0 = torch.tensor([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]).unsqueeze(0)    #(1,3,3)
        factor = torch.tensor([[1., -1., -1.]])     # (1,3)

        c1 = (torch.sum(Gs0 * self.GrPRS, 1)>0).unsqueeze(1)
        c2 = (torch.sum(Gs0 * self.GrPRS, 1)<=0).unsqueeze(1)    # (N,1)
        GsPRS1 = torch.cross(Gr0, self.GrPRS)
        GsPRS2 = torch.cross(Gs0, self.GrPRS)    # (N,3)
        self.GsPRS = c1 * GsPRS1 + c2 * GsPRS2

        GrRad = torch.sum(iR0 * self.GrPRS.unsqueeze(1), 2)
        GsRad = torch.sum(iR0 * self.GsPRS.unsqueeze(1), 2)

        self.GsRad_new = (self.in_plane_rot(GrRad, GsRad)) * factor
        self.GrRad_new = GrRad * factor
        self.GpRad_new = torch.cross(self.GrRad_new, self.GsRad_new)
        
    def in_plane_rot(self, Gr, Gs):
        """Rotate Gs along Gr for angle theta anti-clockwise
        Args:
            Gr:     (N,3)
            Gs:     (N,3)
            theta:  (N,1)
        Returns:
            Gs: Gs after theta rotation, Nx3
        """
        Gp = torch.cross(Gr, Gs)
        Gs_new = torch.cos(self.theta) * Gs + torch.sin(self.theta) * Gp
        return Gs_new

    def rot_baseK(self):
        """Rotate baseK according to Gp,Gr,Gs to generate full k-space trajectory
        Args:
            GpRad:  (N,3)
            GrRad:  (N,3)
            GsRad:  (N,3)
            baseK:  (M,3)
        Returns:
            ktraj:  (3,MxN)
        """
        tmp =   self.GpRad_new.unsqueeze(2) * self.baseKp + \
                self.GrRad_new.unsqueeze(2) * self.baseKr + \
                self.GsRad_new.unsqueeze(2) * self.baseKs       # (N, 3, M)
        self.ktraj = torch.reshape(torch.permute(tmp, (1, 0, 2)), (3, -1))  # beware torch.reshape expanding order starts from last dim


def gen_Grs(Nsegs, Nshots):
    GRCounter = np.arange(0, Nsegs*Nshots)
    GRCounter = np.reshape(GRCounter, (Nshots, Nsegs), order='F')
    GRCounter = np.transpose(GRCounter)
    GRCounter = np.reshape(GRCounter, (Nsegs*Nshots, 1), order='F')
    Phis = [0.465571231876768, 0.682327803828019]
    kz = np.mod(GRCounter * Phis[0], 1) * 2 - 1
    Polar = np.arccos(kz)
    Azi = np.mod(GRCounter * Phis[1], 1) * 2 * np.pi
    Grs = np.concatenate( ( np.sin(Azi)*np.sin(Polar), np.cos(Azi)*np.sin(Polar), np.cos(Polar) ), axis=1)
    return Grs

if __name__ == '__main__':
    Nsegs = 36
    Nshots = 6
    im_size = (64, 64, 64)
    Grs = gen_Grs(Nsegs, Nshots)
    baseK = np.loadtxt("/home/fs0/qijia/scratch/origin_data/cone_dev/phantom/raw_data_2-3-22/ktraj/mjohnson_ktraj_64_200_9980", delimiter='\t')
    
    opt = OptTheta(Grs, baseK, im_size)
    lr = 0.001
    for i in range(100):
        loss = opt()
        loss.backward()
        with torch.no_grad():
            opt.theta -= lr * opt.theta.grad
            opt.theta.grad.zero_()