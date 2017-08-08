import numpy as np
import os

def read_ydata_single(fname,nx,na) :
    """Read data from a single file
    """
    data=np.loadtxt(fname,unpack=True)
    xarr=data[0,::na]
    aarr=data[1,:na]
    warr_2d=(data[2]).reshape([nx,na])
    warr_1d=np.sum((data[3]).reshape([nx,na]),axis=1)/np.sum((data[4]).reshape([nx,na]),axis=1)
    
    return xarr,aarr,warr_2d,warr_1d

def get_run_stats(prefix,nx,na,nsims) :
    """Read data from all files corresponding to the same run and compute their statistics
    """
    dirname=prefix+"_%d"%nx
    fname_stats=dirname+"/wth_stats"
    
    if not os.path.isfile(fname_stats+".npz") :
        x_d,a_d,w_2d_d,w_1d_d=read_ydata_single(dirname+"/wth_voids.txt",nx,na)
        w_2d_ms=[]; w_1d_ms=[]
        for i in np.arange(nsims) :
            print i
            x,a,w2d,w1d=read_ydata_single(dirname+"/wth_mock_%04d.txt"%(i+1),nx,na)
            w_2d_ms.append(w2d)
            w_1d_ms.append(w1d)
        w_2d_ms=np.array(w_2d_ms)
        w_1d_ms=np.array(w_1d_ms)
        w_2d_mean=np.mean(w_2d_ms,axis=0)
        w_2d_std =np.std(w_2d_ms,axis=0)
        w_1d_mean=np.mean(w_1d_ms,axis=0)
        w_1d_std =np.std(w_1d_ms,axis=0)
        w_1d_covar=np.mean(w_1d_ms[:,:,None]*w_1d_ms[:,None,:],axis=0)-w_1d_mean[:,None]*w_1d_mean[:,None]
        np.savez(fname_stats,x=x_d,alpha=a_d,
                 w_2d_data=w_2d_d,w_2d_mean=w_2d_mean,w_2d_error=w_2d_std,
                 w_1d_data=w_1d_d,w_1d_mean=w_1d_mean,w_1d_error=w_1d_std,w_1d_covar=w_1d_covar)
        
    return np.load(fname_stats+".npz")
