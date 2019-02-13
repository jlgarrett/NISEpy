# -*- coding: utf-8 -*-
"""
This code is written to generate the Chlorophyll colormap.
Therefore, instead of importing it, it should just be run.
I'm looking into a more proper way to do this.

Created on Mon Aug 13 2018

@author: joe
"""

import igor.binarywave as igoribw
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

chl  = igoribw.load('Chlorophyll.ibw')
colordata = np.transpose(chl['wave']['wData'])
chlor = {}
L = len(np.transpose(colordata))
chlor['red'] = [(i/(L-1),colordata[0][i]/65533,colordata[0][i]/65533) for i,x in enumerate(colordata[0])]
chlor['green'] = [(i/(L-1),colordata[1][i]/65533,colordata[1][i]/65533) for i,x in enumerate(colordata[0])]
chlor['blue'] = [(i/(L-1),colordata[2][i]/65533,colordata[2][i]/65533) for i,x in enumerate(colordata[0])]
chlorophyll = LinearSegmentedColormap('Chlorophyll', chlor)
plt.register_cmap(cmap=chlorophyll)
