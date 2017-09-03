import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stt
from scipy.interpolate import interp1d
import common as cmm
import sys
from matplotlib.colors import LinearSegmentedColormap

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

#data_milca=cmm.get_run_stats("results_0/output_lens_voids_milca",n_x,n_a,nmocks)
data_milca=cmm.get_run_stats("output_pl60LS_voids_milca",n_x,n_a,nmocks)
data_null =cmm.get_run_stats("output_pl60LS_voids_null" ,n_x,n_a,nmocks)
data_nulln=cmm.get_run_stats("output_pl60LS_voids_nulln",n_x,n_a,nmocks)
data_nilc =cmm.get_run_stats("output_pl60LS_voids_nilc" ,n_x,n_a,nmocks)
data_545  =cmm.get_run_stats("output_pl60LS_voids_545"  ,n_x,n_a,nmocks)
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
    vo=np.ones(len(ids))
    prod_dd=np.dot(vd,np.dot(invcov,vd))
    prod_dt=np.dot(vd,np.dot(invcov,vt))
    prod_do=np.dot(vd,np.dot(invcov,vo))
    prod_tt=np.dot(vt,np.dot(invcov,vt))
    prod_to=np.dot(vt,np.dot(invcov,vo))
    prod_oo=np.dot(vo,np.dot(invcov,vo))
    chi2_null=prod_dd
    chi2_null_sims=np.sum(np.dot(vsims,invcov)*vsims,axis=1)
    a_bf=prod_dt/prod_tt
    of_bf=prod_do/prod_oo
    of_err=np.sqrt(1./prod_oo)
    cov_o=np.linalg.inv(np.array([[prod_tt,prod_to],[prod_to,prod_oo]]))
    a_bfo=np.dot(cov_o,np.array([prod_dt,prod_do]))
    a_err=np.sqrt(1./prod_tt)
    chi2_model=np.dot(vd-a_bf*vt,np.dot(invcov,vd-a_bf*vt))
    chi2_modelo=np.dot(vd-a_bfo[0]*vt-a_bfo[1]*vo,np.dot(invcov,vd-a_bfo[0]*vt-a_bfo[1]*vo))

    d_out={'ndof_null':len(ids),'ndof_model':len(ids)-1, #Number of degrees of freedom
           'of_bf':of_bf,'of_err':of_err, #Best fit amplitude and error
           'a_bf':a_bf,'a_err':a_err, #Best fit amplitude and error
           'chi2_model':chi2_model, #Chi^2 for the model
           'chi2_null':chi2_null, #Chi^2 of the null hypothesis
           'chi2_null_sims':chi2_null_sims, #Null Chi^2 for the mocks
           'p_null_counts':(len(np.where(chi2_null_sims>=chi2_null)[0])+0.)/nmocks, #Null significance from counts
           'p_null':1-stt.chi2.cdf(chi2_null,len(ids)), #Null significance
           'p_model':1-stt.chi2.cdf(chi2_model,len(ids)-1), #Model significance
           'a_bf_o':a_bfo[0],
           'off_bf_o':a_bfo[1],
           'cov_o':cov_o,
           'chi2_model_o':chi2_modelo,
           'p_model_o':1-stt.chi2.cdf(chi2_modelo,len(ids)-2),
           'ndof_model_o':len(ids)-2
           } 

    return d_out

rmax=2.0
st_milca=fit_void(rmax,data_milca)
st_nilc =fit_void(rmax,data_nilc)
st_null =fit_void(rmax,data_null)
st_nulln=fit_void(rmax,data_nulln)
st_545  =fit_void(rmax,data_545)
print " MILCA:"
print " Best fit amplitude : %.3lf += %.3lf"%(st_milca['a_bf'],st_milca['a_err'])
print " Chi^2 = %.3lf, N_dof = %d"%(st_milca['chi2_model'],st_milca['ndof_model'])
print " P(model) : %.3lE"%(st_milca['p_model'])
print " "
print " NILC:"
print " Best fit amplitude : %.3lf += %.3lf"%(st_nilc['a_bf_o'],np.sqrt(st_nilc['cov_o'][0,0]))
print " Best fit offset    : %.3lE += %.3lE"%(st_nilc['off_bf_o'],np.sqrt(st_nilc['cov_o'][1,1]))
print " Chi^2 = %.3lf, N_dof = %d"%(st_nilc['chi2_model_o'],st_nilc['ndof_model_o'])
print " P(model) : %.3lE"%(st_nilc['p_model_o'])
print " "
print " P_MILCA(null hypothesis) : %.3lE"%(st_milca['p_null_counts'])
print " P_NILC(null hypothesis) : %.3lE"%(st_nilc['p_null_counts'])
print " P_NULL(null hypothesis) : %.3lE"%(st_null['p_null_counts'])
print " P_NULLN(null hypothesis) : %.3lE"%(st_nulln['p_null_counts'])
print " P_545(null hypothesis) : %.3lE"%(st_545['p_null_counts'])
print " Chi^2_MILCA(null hypothesis) : %.3lf"%(st_milca['chi2_null'])
print " Chi^2_NILC(null hypothesis) : %.3lf"%(st_nilc['chi2_null'])
print " Chi^2_NULL(null hypothesis) : %.3lf"%(st_null['chi2_null'])
print " Chi^2_NULLn(null hypothesis) : %.3lf"%(st_nulln['chi2_null'])
print " Chi^2_545(null hypothesis) : %.3lf"%(st_545['chi2_null'])
print " "
print " NULL offset (MILCA) : %lE +- %lE"%(st_null['of_bf'],st_null['of_err'])
print " NULL offset (NILC)  : %lE +- %lE"%(st_nulln['of_bf'],st_nulln['of_err'])
print " # sigmas : %.3lf\n"%(np.sqrt(st_milca['chi2_null']-st_milca['chi2_model']))

#Plot covariance matrix
ids=np.where(data_milca['x']<rmax)[0]
covar=((data_milca['w_1d_covar'])[ids,:])[:,ids]
corrmat=covar/np.sqrt((np.diag(covar)[None,:])*(np.diag(covar)[:,None]))
plt.figure()
ax=plt.gca()
im=ax.imshow(corrmat,interpolation='nearest',origin='lower',extent=[0,rmax,0,rmax],
             cmap=plt.get_cmap('bone'),vmin=0,vmax=1)
cb=plt.colorbar(im,ax=ax)
ax.set_xlabel("$\\theta_2/\\theta_v$",fontsize=14)
ax.set_ylabel("$\\theta_1/\\theta_v$",fontsize=14)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(12)
plt.savefig("doc/corrmat.pdf",bbox_inches='tight')

#Plot fiducial results
plt.figure()
ax=plt.gca()
ax.errorbar(data_milca['x'],data_milca['w_1d_data']-data_milca['w_1d_mean'],
             yerr=data_milca['w_1d_error'],fmt='ro',lw=2,elinewidth=2,label='$y_{\\rm data}$')
ax.plot(x_th_f,w_th_f                 ,'k--',lw=2,label='$y_{\\rm model}$')
ax.plot(x_th_f,w_th_f*st_milca['a_bf'],'k-',lw=2,label='$\\alpha_v\\,y_{\\rm model}$')
ax.plot([0,2],[0,0],'k--',lw=1)
ax.text(0.1,4.5E-8,'$\\alpha_v=%.2lf\\pm%.2lf$'%(st_milca['a_bf'],st_milca['a_err']),fontsize=18)
ax.set_ylim([-6.2E-8,6.2E-8])
ax.set_xlim([0,2])
ax.set_xlabel("$\\theta/\\theta_v$",fontsize=14)
ax.set_ylabel("$\\bar{y}(\\theta)$",fontsize=14)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(12)
plt.legend(loc='lower right',frameon=False,fontsize=16)#,ncol=2)
plt.savefig("doc/y_result.pdf",bbox_inches='tight')

#Systematics plot
plt.figure()
ax=plt.gca()
y1=np.maximum((data_545['w_1d_data']-data_545['w_1d_mean'])*(alpha_cib+s_alpha_cib),(data_545['w_1d_data']-data_545['w_1d_mean'])*(alpha_cib-s_alpha_cib))
y2=np.minimum((data_545['w_1d_data']-data_545['w_1d_mean'])*(alpha_cib+s_alpha_cib),(data_545['w_1d_data']-data_545['w_1d_mean'])*(alpha_cib-s_alpha_cib))
xext=np.linspace(data_545['x'][0],data_545['x'][-1],len(data_545['x'])-1)
y1f=interp1d(data_545['x'],y1); y1ext=y1f(xext)
y2f=interp1d(data_545['x'],y2); y2ext=y2f(xext)
ax.fill_between(xext,y1ext,y2ext,facecolor='#AAAAAA')
ax.errorbar(data_545['x'],(data_545['w_1d_data']-data_545['w_1d_mean'])*alpha_cib,yerr=data_545['w_1d_error']*alpha_cib,
            fmt='kD',lw=2,elinewidth=2,ms=4,label='${\\rm Leakage}$')
ax.errorbar(data_milca['x'],data_milca['w_1d_data']-data_milca['w_1d_mean'],
             yerr=data_milca['w_1d_error'],fmt='ro',lw=2,elinewidth=2,label='$y_{\\rm MILCA}$')
ax.errorbar(data_nilc['x']+0.02,data_nilc['w_1d_data']-data_nilc['w_1d_mean']-st_nilc['off_bf_o'],
             yerr=data_nilc['w_1d_error'],fmt='bs',lw=2,elinewidth=2,label='$y_{\\rm NILC}$')
ax.plot(x_th_f,w_th_f*st_milca['a_bf'],'k-',lw=2,label='${\\rm best\\,\\,fit\\,\\,MILCA}$')
ax.plot([0,2],[0,0],'k--',lw=1)
ax.set_ylim([-6.2E-8,6.2E-8])
ax.set_xlim([0,2])
ax.set_xlabel("$\\theta/\\theta_v$",fontsize=14)
ax.set_ylabel("$\\bar{y}(\\theta)$",fontsize=14)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(12)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(12)
plt.legend(loc='lower right',frameon=False,fontsize=14,ncol=1)
plt.savefig("doc/y_syst.pdf",bbox_inches='tight')

#Plot 2D stacks
np.random.seed(1234)
a_refac=8
aarr=np.radians(np.linspace(0,360,n_a*a_refac))
xarr=data_milca['x']
r,theta=np.meshgrid(xarr,aarr)

cdict={'red':  ((0.00,0.00,0.00),
                (0.25,0.20,0.20),
                (0.50,1.00,1.00),
                (0.75,1.00,1.00),
                (1.00,0.50,0.50)),
       'green':((0.00,0.00,0.00),
                (0.25,0.20,0.20),
                (0.50,1.00,1.00),
                (0.75,0.40,0.40),
                (1.00,0.20,0.20)),
       'blue' :((0.00,0.50,0.50),
                (0.25,1.00,1.00),
                (0.50,1.00,1.00),
                (0.75,0.00,0.00),
                (1.00,0.00,0.00))}
cm=LinearSegmentedColormap('BlueRed', cdict)

def plot2d(arr,vmn,vmx) :
    fig,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2,subplot_kw=dict(projection='polar'),figsize=(10,10));

    ipl=0
    for arr,ax in zip(arr,[ax1,ax2,ax3,ax4]) :
        print ipl
        arrplot=np.zeros([n_a*a_refac,n_x])
        for i in np.arange(n_a) :
            arrplot[i*a_refac:(i+1)*a_refac,:]=(arr[:,i])[None,:]
        p=ax.pcolor(theta,r,arrplot,vmin=vmn,vmax=vmx,cmap=cm)#plt.get_cmap('seismic'));
        ipl+=1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(p, cax=cbar_ax)
    plt.savefig("doc/stacks_2d.pdf",bbox_inches='tight')

plot2d([data_milca['w_2d_data']-data_milca['w_2d_mean'],
        data_milca['w_2d_sims'][np.random.randint(nmocks)]-data_milca['w_2d_mean'],
        data_milca['w_2d_sims'][np.random.randint(nmocks)]-data_milca['w_2d_mean'],
        data_milca['w_2d_sims'][np.random.randint(nmocks)]-data_milca['w_2d_mean']],
       -6.2E-8,6.2E-8)
plt.show()
