# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 12:52:58 2022

@author: Chengjie Luo
"""



from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

outputfolder=r"C:\Users\s_nab\OneDrive\Documents\BEP\PhysiCell\run1"
flag_series=1
tmax=125
dt=1

xmin=-200
xmax=200
ymin=-200
ymax=200
zmin=-200
zmax=200

if flag_series==1:
    first_index=1
    last_index=int(tmax/dt)
    nindex=last_index-first_index+1
    
    
    xs=np.zeros(nindex)
    ys=np.zeros(nindex)
    zs=np.zeros(nindex)
    Nc=np.zeros(nindex)
    
    
    
    
    times=np.zeros(nindex)
    
    f = open(outputfolder+"/cells_series.xyz", "w")
    fd=open(outputfolder+"/cells_series.dump","w")
    for i,n in enumerate(range( first_index,last_index+1,1)):
        try:
            filename='output'+"%08i"%n+'.xml'
            print(filename)
            mcds=pyMCDS(filename,outputfolder)
            times[i]= mcds.get_time()
            Nc[i]=mcds.data['discrete_cells']['position_x'].shape[0]
            
            ID=mcds.data['discrete_cells']['ID']
            x=mcds.data['discrete_cells']['position_x']
            y=mcds.data['discrete_cells']['position_y']
            z=mcds.data['discrete_cells']['position_z']
            v=mcds.data['discrete_cells']['total_volume']
            celltype=mcds.data['discrete_cells']['cell_type']
            


            

            ###### write to dump format
            fd.write("ITEM: TIMESTEP \n")
            fd.write(str(int(times[i]))+'\n')
            fd.write("ITEM: NUMBER OF ATOMS \n")
            fd.write(str(int(Nc[i]))+'\n')
            fd.write("ITEM: BOX BOUNDS pp pp pp\n")
            fd.write(str(xmin)+' '+str(xmax)+ "\n")
            fd.write(str(ymin)+' '+str(ymax)+ "\n")
            fd.write(str(zmin)+' '+str(zmax)+ "\n")
            fd.write("ITEM: ATOMS type id x y z radius \n")
            for ic in np.arange(int(Nc[i])):
                ###### write to xyz
                tmpx=str(x[ic])
                tmpy=str(y[ic])
                tmpz=str(z[ic])
                tmpid=str(int(ID[ic]))
                tmpv=v[ic]
                tmpr=str((tmpv/(4/3.*np.pi))**(1/3))
                tmptype=str(celltype[ic]+1)
                
                fd.write(tmptype+' '+tmpid+' '+tmpx+' '+tmpy+' '+tmpz+' '+tmpr+'\n')
        except:
            print("not all included")     
            break
    
    
    
    
    f.close()
    fd.close()
    
    print(Nc)
    


