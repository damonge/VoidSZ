import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stt
from scipy.interpolate import interp1d
import common as cmm
import sys

if len(sys.argv)!=2 :
    print "Usage: analysis.py n_x"
    exit(1)
n_x=int(sys.argv[1])
n_a=16
nmocks=1000

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

data_milca=cmm.get_run_stats("output_voids_milca",n_x,n_a,nmocks)
data_null =cmm.get_run_stats("output_voids_null" ,n_x,n_a,nmocks)
data_nilc =cmm.get_run_stats("output_voids_nilc" ,n_x,n_a,nmocks)
data_545  =cmm.get_run_stats("output_voids_545"  ,n_x,n_a,nmocks)
x_th_f,w_th_f,dum=np.loadtxt("data/data_y/y_th_void.txt",unpack=True)
wthf=interp1d(x_th_f,w_th_f)

alpha_cib,s_alpha_cib,alpha_gal,s_alpha_gal=np.loadtxt("data/data_y/y_contamination.txt")
w_th=wthf(data_milca['x'])

def fit_void(xmx,data) :
    ids=np.where(data['x']<=xmx)[0]
    covar=((data['w_1d_covar'])[ids,:])[:,ids]
    invcov=np.linalg.inv(covar)*(nmocks-len(ids)-2.)/(nmocks-1.) #Inverse covariance with correctin factor from Hartlap et al. 2007
    vd=(data['w_1d_data']-data['w_1d_mean'])[ids]
    vsims=(data['w_1d_sims']-data['w_1d_mean'][None,:])[:,ids]
    vt=(wthf(data['x']))[ids]
    prod_dd=np.dot(vd,np.dot(invcov,vd))
    prod_dt=np.dot(vd,np.dot(invcov,vt))
    prod_tt=np.dot(vt,np.dot(invcov,vt))
    chi2_null=prod_dd
    chi2_null_sims=np.sum(np.dot(vsims,invcov)*vsims,axis=1)
    a_bf=prod_dt/prod_tt
    a_err=np.sqrt(1./prod_tt)
    chi2_model=np.dot(vd-a_bf*vt,np.dot(invcov,vd-a_bf*vt))

    d_out={'ndof_null':len(ids),'ndof_model':len(ids)-1, #Number of degrees of freedom
           'a_bf':a_bf,'a_err':a_err, #Best fit amplitude and error
           'chi2_model':chi2_model, #Chi^2 for the model
           'chi2_null':chi2_null, #Chi^2 of the null hypothesis
           'chi2_null_sims':chi2_null_sims, #Null Chi^2 for the mocks
           'p_null_counts':(len(np.where(chi2_null_sims>=chi2_null)[0])+0.)/nmocks, #Null significance from counts
           'p_null':1-stt.chi2.cdf(chi2_null,len(ids)), #Null significance
           'p_model':1-stt.chi2.cdf(chi2_model,len(ids)-1) #Model significance
           } 

    return d_out

st_milca=fit_void(2.,data_milca)
st_null =fit_void(2.,data_null)
print " Best fit amplitude : %.3lf += %.3lf"%(st_milca['a_bf'],st_milca['a_err'])
print " Chi^2 = %.3lf, N_dof = %d"%(st_milca['chi2_model'],st_milca['ndof_model'])
print " P(model) : %.3lE"%(st_milca['p_model'])
print " P(null hypothesis) : %.3lE"%(st_milca['p_null'])

plt.figure()
ax=plt.gca()
ax.plot(data_milca['x'],data_milca['w_1d_sims'][0]-data_milca['w_1d_mean'],'-',color='b',lw=2,alpha=1,label='$y_{\\rm mocks}$')
for w in data_milca['w_1d_sims'] : 
    ax.plot(data_milca['x'],w-data_milca['w_1d_mean'],'-',color='b',lw=2,alpha=0.1)
ax.errorbar(data_milca['x'],data_milca['w_1d_data']-data_milca['w_1d_mean'],
             yerr=data_milca['w_1d_error'],fmt='r-',lw=2,elinewidth=2,label='$y_{\\rm data}$')
ax.plot(x_th_f,w_th_f*st_milca['a_bf'],'k-',lw=2,label='$y_{\\rm theory},\\,\\,{\\rm scaled}$')
ax.plot(x_th_f,w_th_f                 ,'k--',lw=2)#,label='y_{\\rm theory},\\,\\,{\\rm scaled}')
ax.plot(data_545['x'],(data_545['w_1d_data']-data_545['w_1d_mean'])*s_alpha_cib,'y-',lw=2,label='${\\rm CIB\\,\\,leakage}$')
ax.set_ylim([-7.5E-8,7.5E-8])
ax.set_xlim([data_milca['x'][0],2])
ax.set_xlabel("$\\theta/\\theta_v$",fontsize=14)
ax.set_ylabel("$\\langle y(\\theta)\\rangle$",fontsize=14)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(12)
plt.legend(loc='lower right',frameon=False,fontsize=14)
plt.savefig("doc/y_result.pdf",bbox_inches='tight')
plt.show()

exit(1)

aarr=np.radians(np.linspace(0,360,n_a))
xarr=x_m
r,theta=np.meshgrid(xarr,aarr)

def plot2d(arr,vmn,vmx) :
    fig,ax=plt.subplots(subplot_kw=dict(projection='polar'));
    p=ax.pcolor(theta,r,np.transpose(arr),vmin=vmn,vmax=vmx);
    cb=plt.colorbar(p,ax=ax)
plot2d(w_2d_d      ,-1E-7,5E-8)
plot2d(w_2d_n      ,-1E-7,5E-8)
plot2d(w_2d_ms[ 0] ,-1E-7,5E-8)
plot2d(w_2d_m      ,-1E-7,5E-8)
plot2d(ew_2d_m     , 0   ,6E-8)
