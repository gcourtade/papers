#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np


# In[6]:


csi3  = np.loadtxt('AA14_peptide_assignment.out', skiprows=1, dtype=str)
csi3


# In[10]:


resnr = csi3[:,0].astype(int)
rci = np.where(csi3[:,4] == 'NA','0.0', csi3[:,4]).astype(float)


# In[17]:


# A lower RCI value indicates a more rigid region
# This paper should be additionally cited, e.g. in the figure caption, when using RCI from CSI3.0:
# https://pubs.acs.org/doi/10.1021/ja054842f
#
plt.plot(resnr, rci, 'k-', label='RCI from CSI3')
plt.ylabel('RCI from CSI3.0')
plt.xlabel('Residue Number')
plt.xlim(resnr.min(), resnr.max())


# In[ ]:




