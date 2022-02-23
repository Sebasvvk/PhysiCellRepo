# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 20:23:08 2022

@author: njulu
"""



from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

plotfromtxt=1
adhesions=[0,0.4,0.8,1.6]
speeds=[0,0.05,0.1,0.3,0.5,0.8]
cycles=[0.00018,0.00036,0.00072,0.00144,0.00288,0.005,0.01]
#follower_speeds=[0.1,0.2,0.3,0.4]
for ia, cycle in enumerate(cycles):
    for isp, adhesion in enumerate(adhesions):
        
        plt.figure(ia*100+isp*1000+1)
        plt.loglog(np.arange(1,100),489*np.exp(0.00036*np.arange(1,100)*60),'k--')
        plt.loglog(np.arange(1,100),489*np.exp(0.00072*np.arange(1,100)*60),'k--')
        
        plt.figure(ia*100+isp*1000+2)
        plt.loglog(np.arange(1,100),410*np.arange(1,100),'k--')
        
        


        for speed in speeds:
            pre='onesphere_adh_'+str(adhesion)+'_speed_'+str(speed)+'_cycle_'+str(cycle)
            folder_name='./output_'+pre
            
            
            if (plotfromtxt==1):
                data=np.loadtxt(folder_name+'/N_t.txt')
                plt.figure(ia*100+isp*1000+1)
                plt.plot(data[:,0],data[:,1],'o',label='adh_'+str(adhesion)+'_speed_'+str(speed)+'_cycle_'+str(cycle))
                data=np.loadtxt(folder_name+'/msd.txt')
                plt.figure(ia*100+isp*1000+2)
                plt.plot(data[:,0],data[:,1],label='adh_'+str(adhesion)+'_speed_'+str(speed)+'_cycle_'+str(cycle))
                        
                
            else:
            
            
            
            
                outputfolder=folder_name
                flag_series=1
                tmax=125
                dt=1
                
                xmin=-300
                xmax=300
                ymin=-300
                ymax=300
                zmin=-300
                zmax=300
                
                if flag_series==1:
                    first_index=1
                    last_index=int(tmax/dt)
                    nindex=last_index-first_index+1
                    
                    
                    xs=np.zeros(nindex)
                    ys=np.zeros(nindex)
                    zs=np.zeros(nindex)
                    Nc=np.zeros(nindex)
                    
                    
                    
                    
                    times=np.zeros(nindex)
                    ts=[]
                    Ns=[]
                    msds=[]
                
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
                            if (i==0):
                                x0=x
                                y0=y
                                z0=z
                                ID0=ID
                                v0=v
                                celltype0=celltype
                            
                            else:
                                msd=np.mean((x[:x0.shape[0]]-x0)**2+(y[:x0.shape[0]]-y0)**2+(z[:x0.shape[0]]-z0)**2)
                                msds.append(msd)
                                ts.append(i)
                                Ns.append(ID.shape[0])
                                
                
                
                            
                
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
                    
                    ts=np.array(ts)
                    Ns=np.array(Ns)
                    msds=np.array(msds)
                    data=np.concatenate((ts.reshape(-1,1),Ns.reshape(-1,1)),axis=1)
                    np.savetxt(folder_name+'/N_t.txt',data)
                    data=np.concatenate((ts.reshape(-1,1),msds.reshape(-1,1)),axis=1)
                    np.savetxt(folder_name+'/msd.txt',data)
                    plt.figure(ia*100+isp*1000+1)
                    plt.plot(ts,Ns,'o',label='adh_'+str(adhesion)+'_speed_'+str(speed)+'_cycle_'+str(cycle))
                    plt.legend()
                    plt.ylabel('N')
                    plt.xlabel('t')
                    plt.figure(ia*100+isp*1000+2)
                    plt.plot(ts,msds,'o-',label='adh_'+str(adhesion)+'_speed_'+str(speed)+'_cycle_'+str(cycle))
                    plt.legend()
                    plt.ylabel('msd')
                    plt.xlabel('t')
                    
                
                fd.close()
            plt.legend()
        plt.show()
    


