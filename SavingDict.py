# -*- coding: utf-8 -*-
"""
Created on Sat Mar 24 13:17:40 2018

@author: joe
"""

import tables
import CypherKPFM as ck

class Coordinates(tables.IsDescription):
    x_coordinate = tables.Int32Col()
    y_coordinate = tables.Int32Col()
    
    #testfile.create_array(group2,dict24nm[k[0]]['cleanlabels']

def createTextDict( original_dict ):
    keys = list(original_dict.keys())
    keys.sort
    new_dict = {}
    for i in keys:
        if len(original_dict[i]) < 16:
            new_dict[i] = tables.StringCol(itemsize=16)
    
    return new_dict
    
class TanhFits(tables.IsDescription):
    x_coordinate = tables.Int32Col()
    y_coordinate = tables.Int32Col()
    Vb = tables.Float32Col()
    a = tables.Float32Col()
    Voffset = tables.Float32Col()
    x0 = tables.Float32Col()
    chisq = tables.Float32Col()
    Vb_err = tables.Float32Col()
    a_err = tables.Float32Col()
    Voffset_err = tables.Float32Col()
    x0_err = tables.Float32Col()
    
def loadonefit( onefit, KPFMdictionary):
    myKeys = list(KPFMdictionary['params'].keys())
    myKeys.sort()
    L = len(myKeys)
    
    for i in range(L):
        onefit['x_coordinate'] = myKeys[i][0]
        onefit['y_coordinate'] = myKeys[i][1]
        onefit['Vb'] = KPFMdictionary['params'][myKeys[i]][0][0]
        onefit['a'] = KPFMdictionary['params'][myKeys[i]][0][1]
        onefit['Voffset'] = KPFMdictionary['params'][myKeys[i]][0][2]
        onefit['x0'] = KPFMdictionary['params'][myKeys[i]][0][3]
        onefit['chisq'] = KPFMdictionary['params'][myKeys[i]][1]
        onefit['Vb_err'] = KPFMdictionary['uncertainty'][myKeys[i]][0]
        onefit['a_err'] = KPFMdictionary['uncertainty'][myKeys[i]][1]
        onefit['Voffset_err'] = KPFMdictionary['uncertainty'][myKeys[i]][2]
        onefit['x0_err'] = KPFMdictionary['uncertainty'][myKeys[i]][3]
        onefit.append()
    
    #Remember to flush!
    return 0

def saveWholeDict( d, name ):
    with tables.open_file(name+'.h5', mode = 'w', title = name) as savingFile:
        k = list(d.keys())
        k.sort()
        for i in k:
            group = savingFile.create_group('/', i, i + ' data')
            imagegroup = savingFile.create_group('/'+i,
                                             'imageData', 'Saved Images')
            labels = list(d[i]['data']['cleanlabels'].keys())
            labels.sort()
            for j in labels:  
                savingFile.create_array(imagegroup,j,
                                  ck.get_image(d[i]['data'],j),j) 
        
            newdict = createTextDict(d[i]['data']['cleannote'])
            table2 = savingFile.create_table(group, "note",
                                         newdict, "Metadata note")
        
            ndk = list(newdict.keys())
            ndk.sort()
            note_attributes = table2.row
            for j in ndk:
                note_attributes[j] = d[i]['data']['cleannote'][j]
        
            note_attributes.append()
            table2.flush()
        
            fitstable = savingFile.create_table(group, "FitParameters",
                                                TanhFits, "Fitting Parameters")
        
            onefit = fitstable.row
            loadonefit(onefit, d[i])
            fitstable.flush()
            gcuts = savingFile.create_group(group,
                                      "cuts", "Voltage versus distance cuts")
   
            myKeys = list(d[i]['cuts'].keys())
            myKeys.sort()
            for j in range(len(myKeys)):
                savingFile.create_array(gcuts, "cut_{}".format(j),
                        d[i]['cuts'][myKeys[j]],
                        'cut centered on pixel {}'.format(str(myKeys[j])))
        
        savingFile.close()
    return 0
        