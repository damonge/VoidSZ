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

data_milca=cmm.get_run_stats("output_voids_milca",n_x,n_a,nmocks)
data_null =cmm.get_run_stats("output_voids_null" ,n_x,n_a,nmocks)
data_nilc =cmm.get_run_stats("output_voids_nilc" ,n_x,n_a,nmocks)
data_545  =cmm.get_run_stats("output_voids_545"  ,n_x,n_a,nmocks)
alpha_cib,s_alpha_cib,alpha_gal,s_alpha_gal=np.loadtxt("data/data_y/y_contamination.txt")

plt.errorbar(data_milca['x'],data_milca['w_1d_data']-data_milca['w_1d_mean'],
             yerr=data_milca['w_1d_error'],color='r')
plt.errorbar(data_null['x'],data_null['w_1d_data']-data_null['w_1d_mean'],
             yerr=data_null['w_1d_error'],color='y')
plt.errorbar(data_nilc['x'],data_nilc['w_1d_data']-data_nilc['w_1d_mean'],
             yerr=data_nilc['w_1d_error'],color='b')
plt.errorbar(data_545['x'],(data_545['w_1d_data']-data_545['w_1d_mean'])*s_alpha_cib,
             yerr=data_545['w_1d_error']*s_alpha_cib,color='k')
plt.show()
exit(1)

def read_sim(fname,nx,na) :
    data=np.loadtxt(fname,unpack=True)
    xarr=data[0,::na]
    aarr=data[0,:na]
    warr_2d=(data[2]).reshape([nx,na])
    warr_1d=np.sum((data[3]).reshape([nx,na]),axis=1)/np.sum((data[4]).reshape([nx,na]),axis=1)
    
    return xarr,aarr,warr_2d,warr_1d

x_d,a_d,w_2d_d,w_1d_d=read_sim("output_voids_%d/wth_voids.txt"%n_x,n_x,n_a)
x_n,a_n,w_2d_n,w_1d_n=read_sim("output_voids_null_%d/wth_voids.txt"%n_x,n_x,n_a)
#x_c,a_c,w_2d_c,w_1d_c=read_sim("output_voids_545_%d/wth_voids.txt"%n_x,n_x,n_a)
x_ms=[]; a_ms=[]; w_2d_ms=[]; w_1d_ms=[]
#x_msc=[]; a_msc=[]; w_2d_msc=[]; w_1d_msc=[]
for i in np.arange(nmocks) :
    x,a,w2d,w1d=read_sim("output_voids_%d/"%n_x+"wth_mock_%04d.txt"%(i+1),n_x,n_a)
    x_ms.append(x)
    a_ms.append(a)
    w_2d_ms.append(w2d)
    w_1d_ms.append(w1d)
#    x,a,w2d,w1d=read_sim("output_voids_545_%d/"%n_x+"wth_mock_%04d.txt"%(i+1),n_x,n_a)
#    x_msc.append(x)
#    a_msc.append(a)
#    w_2d_msc.append(w2d)
#    w_1d_msc.append(w1d)
x_ms=np.array(x_ms)
a_ms=np.array(a_ms)
w_2d_ms=np.array(w_2d_ms)
w_1d_ms=np.array(w_1d_ms)
x_m=np.mean(x_ms,axis=0)
a_m=np.mean(a_ms,axis=0)
w_2d_m=np.mean(w_2d_ms,axis=0)
w_1d_m=np.mean(w_1d_ms,axis=0)
ew_2d_m=np.std(w_2d_ms,axis=0)
ew_1d_m=np.std(w_1d_ms,axis=0)
#x_msc=np.array(x_msc)
#a_msc=np.array(a_msc)
#w_2d_msc=np.array(w_2d_msc)
#w_1d_msc=np.array(w_1d_msc)
#x_mc=np.mean(x_msc,axis=0)
#a_mc=np.mean(a_msc,axis=0)
#w_2d_mc=np.mean(w_2d_msc,axis=0)
#w_1d_mc=np.mean(w_1d_msc,axis=0)
#ew_2d_mc=np.std(w_2d_msc,axis=0)
#ew_1d_mc=np.std(w_1d_msc,axis=0)

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
#plot2d(w_2d_mc     ,None,None)
#plot2d(ew_2d_mc    ,None,None)

#plt.figure(); plt.errorbar(x_c,w_1d_c-w_1d_mc,yerr=ew_1d_mc);

x_th_f,w_th_f,xi_th_f=np.loadtxt("y_th_void.txt",unpack=True)
wthf=interp1d(x_th_f,w_th_f)
w_th=wthf(x_d)

def get_chis(xmx,plot_stuff=False) :
    ids=np.where(x_d<=xmx)[0]
    covar=np.mean(w_1d_ms[:,ids,None]*w_1d_ms[:,None,ids],axis=0)-w_1d_m[ids,None]*w_1d_m[None,ids]
    if plot_stuff :
        plt.figure();
        plt.imshow(covar/np.sqrt(np.diag(covar)[:,None]*np.diag(covar)[None,:]),interpolation='nearest',origin='lower'); 
    invcov=np.linalg.inv(covar)
    vd=w_1d_d-w_1d_m
    vn=w_1d_n-w_1d_m
    vt=w_th
    chi2d=np.dot(vd[ids],np.dot(invcov,vd[ids]))
    chi2n=np.dot(vn[ids],np.dot(invcov,vn[ids]))
    chi2ms=np.array([np.dot((w-w_1d_m)[ids],np.dot(invcov,(w-w_1d_m)[ids])) for w in w_1d_ms])
    pd=(len(np.where(chi2ms>=chi2d)[0])+0.)/nmocks
    pn=(len(np.where(chi2ms>=chi2n)[0])+0.)/nmocks
    ic_tt=np.dot(vt[ids],np.dot(invcov,vt[ids]))
    ic_dt=np.dot(vd[ids],np.dot(invcov,vt[ids]))
    bbf=1.
    chi2t=np.dot((vd-bbf*vt)[ids],np.dot(invcov,(vd-bbf*vt)[ids]))
    bbf=ic_dt/ic_tt
    berr=1./np.sqrt(ic_tt)
    pt=(len(np.where(chi2ms>=chi2t)[0])+0.)/nmocks

    return chi2d,chi2n,chi2t,chi2ms,pd,pn,pt,bbf,berr,len(ids)

chi2_d,chi2_n,chi2_t,chi2_ms,prob,probn,probt,b_bf,b_err,ndof=get_chis(2.,plot_stuff=True)

plt.figure()
plt.plot(x_m,w_1d_ms[0]-w_1d_m,'-',color='b',alpha=1,lw=2,label='$y_{\\rm mocks}$')
for w in w_1d_ms :
    plt.plot(x_m,w-w_1d_m,'-',color='b',alpha=0.1,lw=2)
plt.plot(x_d,w_1d_n-0*w_1d_m,'y-',lw=2,label='$y_{\\rm null}$')
plt.errorbar(x_d,w_1d_d-w_1d_m,yerr=ew_1d_m*np.sqrt(920./774),
             fmt='r-',lw=2,elinewidth=2,label='$y_{\\rm data}$')
plt.plot(x_th_f,b_bf*w_th_f,'k-',lw=2,label='$y_{\\rm theory}$')
plt.plot(x_th_f,w_th_f,'k--')#,lw=2,label='$y_{\\rm theory}$')
plt.xlabel('$r/r_{\\rm void}$',fontsize=18)
plt.ylabel('$\\langle y(r)\\rangle$',fontsize=18)
plt.ylim([-7.5E-8,7.5E-8])
plt.xlim([0,2])
plt.legend(loc='lower right',frameon=False,fontsize=18)
plt.savefig('y_voids.png',bbox_inches='tight')

fracs=0.2*(np.arange(15)+1)
ndofs=np.zeros_like(fracs)
chi2s=np.zeros_like(fracs)
probs=np.zeros_like(fracs)
probsn=np.zeros_like(fracs)
probsb=np.zeros_like(fracs)
for i in np.arange(len(fracs)) :
    cd,cn,ct,cm,pr,prn,prt,bb,be,nd=get_chis(fracs[i])
    ndofs[i]=nd+0.
    chi2s[i]=cd
    probs[i]=pr
    probsn[i]=prn
    probsb[i]=1-stt.chi2.cdf(cd,nd)

print "chi^2(0) = %.3lf, "%chi2_d+"chi^2(null) = %.3lf, "%chi2_n+"<chi^2(mocks)> = %.3lf"%(np.mean(chi2_ms))
print "b_BF = %.3lf +-"%b_bf+" %.3lf"%b_err
print "S/N = %.3lf"%(b_bf/b_err)
print "chi^2 = %.3lf, "%chi2_t+"n_dof = %d, "%(ndof-1)+"PTE = %.3lE"%(1-stt.chi2.cdf(chi2_t,ndof-1))
#plt.figure(); plt.plot(fracs,chi2s/ndofs)
#plt.figure(); plt.plot(fracs,probs,'r-'); plt.plot(fracs,probsb,'g-'); plt.plot(fracs,probsn,'b-')
#plt.figure(); plt.hist(chi2_ms,bins=50)
plt.show()
