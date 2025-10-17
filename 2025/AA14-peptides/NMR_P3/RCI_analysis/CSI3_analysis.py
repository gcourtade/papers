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


# In[49]:


# A lower RCI value indicates a more rigid region
# This paper should be additionally cited, e.g. in the figure caption, when using RCI from CSI3.0:
# https://pubs.acs.org/doi/10.1021/ja054842f
#
plt.plot(resnr, rci, 'k.-', markersize=5, linewidth=1.6)
plt.plot(resnr[7:14], rci[7:14], '.-', markersize=8, linewidth=2, color='#f1bc3d', label='Helical region')
plt.ylabel('RCI from CSI3.0')
plt.xlabel('Residue Number')
plt.xlim(resnr.min(), resnr.max())
plt.legend()
plt.savefig('AA14_peptide_RCI_CSI3.png', dpi=300)
plt.savefig('AA14_peptide_RCI_CSI3.pdf', dpi=300)
plt.show()


# In[ ]:




