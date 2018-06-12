# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 12:01:00 2018

@author: joe
"""

def Make_test_KPFM_data():
    KPFMFile = {}
    
    KPFMFile['cleanlabels'] = {'test': 0}
    L=8
    KPFMFile['cleannote'] = {'dx': '1', 'Initial PointsLines': str(L), 'Initial ScanPoints':str(L)}
    practicedata = np.zeros((L,L,0))
    #KPFMFile['cleanlabels'] = {'UserIn1Retrace':0, 'UserIn2Retrace':2,'UserIn1Trace':1,'UserIn2Trace':3}
    KPFMFile['wave'] = {'wData':practicedata}
    
    ui1t = np.zeros((L,L))-.1
    ui1r = np.zeros((L,L))-.1
    ui2t = np.zeros((L,L))
    ui2r = np.zeros((L,L))
    
    ui2t = ui2t - 1.8e-4
    ui2r = ui2r - 1.8e-4
    ui2r[1,1] = 3e-3
    ui1r[1,1] = .1
    
    ui2r[2:6,3:6] =3e-3
    ui2t[4,5:7] =-3e-3
    ui1t[4,5:7] = .05
    ui1r[1:6,3:6] = -.15
    
    
    #max(newfile['cleanlabels'].values()) + 1
    add_image( KPFMFile, ui2t ,'UserIn2Trace')
    add_image( KPFMFile, ui2r ,'UserIn2Retrace')
    add_image( KPFMFile, ui1t ,'UserIn1Trace')
    add_image( KPFMFile, ui1r ,'UserIn1Retrace')
    
    return KPFMFile