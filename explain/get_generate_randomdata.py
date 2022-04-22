#!/usr/bin/env python
import sys
import numpy as np
prefix=sys.argv[1]
win_length=128000
dim=4
mask=np.ones((win_length,dim))
template=np.zeros((win_length,dim))
y_mask_rela_start=450
y_mask_rela_end=550
mask_rela_start=y_mask_rela_start*128
mask_rela_end=y_mask_rela_end*128
y=np.zeros(1000)
y[y_mask_rela_start:y_mask_rela_end]=1000000
np.savez_compressed("%s.npz"%(prefix),mask=mask,template=template,y_scale=y,y_real=np.zeros(1000),y_mask_start=y_mask_rela_start,y_mask_end=y_mask_rela_end,X_mask_start=mask_rela_start,X_mask_end=mask_rela_end)