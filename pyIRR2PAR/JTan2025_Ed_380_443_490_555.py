import numpy as np
from copy import deepcopy
import os

def b_spline_basis(x):
    """
    tool to generate b-spline basis using vectorized De Boor recursion
    the basis functions extrapolate linearly past the end-knots.

    code is based on:
    Servén D., Brummitt C. (2018). pyGAM: Generalized Additive Models in Python. Zenodo. 
    DOI: 10.5281/zenodo.1208723
    
    Parameters
    ----------
    x : array-like, with ndims == 1.

    Returns
    -------
    basis : array containing b-spline basis functions
            with shape (len(x), n_splines)
    """

    edge_knots =  np.array([0,200.])
    n_splines = 100
    spline_order = 3
    if np.ravel(x).ndim != 1:
        raise ValueError('Data must be 1-D, but found {}'.format(np.ravel(x).ndim))

    # rescale edge_knots to [0,1], and generate boundary knots
    edge_knots = np.sort(deepcopy(edge_knots))
    offset = edge_knots[0]
    scale = edge_knots[-1] - edge_knots[0]
    
    boundary_knots = np.linspace(0, 1, 1 + n_splines - spline_order)
    diff = np.diff(boundary_knots[:2])[0]
    
    # rescale x as well
    x = (np.ravel(deepcopy(x)) - offset) / scale

    # append 0 and 1 in order to get derivatives for extrapolation
    x = np.r_[x, 0.0, 1.0]

    # determine extrapolation indices
    x_extrapolte_l = x < 0
    x_extrapolte_r = x > 1
    x_interpolate = ~(x_extrapolte_r + x_extrapolte_l)

    # formatting
    x = np.atleast_2d(x).T

    # augment knots
    aug = np.arange(1, spline_order + 1) * diff
    aug_knots = np.r_[-aug[::-1], boundary_knots, 1 + aug]
    aug_knots[-1] += 1e-9  # want last knot inclusive

    # prepare Haar Basis
    bases = (x >= aug_knots[:-1]).astype(int) * (x < aug_knots[1:]).astype(int)
    bases[-1] = bases[-2][::-1]  # force symmetric bases at 0 and 1

    # do recursion from Hastie et al. vectorized
    maxi = len(aug_knots) - 1
    for m in range(2, spline_order + 2):
        maxi -= 1

        # left sub-basis
        num = x - aug_knots[:maxi]
        num *= bases[:, :maxi]
        denom = aug_knots[m - 1 : maxi + m - 1] - aug_knots[:maxi]
        left = num / denom

        # right sub-basis
        num = (aug_knots[m : maxi + m] - x) * bases[:, 1 : maxi + 1]
        denom = aug_knots[m : maxi + m] - aug_knots[1 : maxi + 1]
        right = num / denom

        # track previous bases and update
        prev_bases = bases[-2:]
        bases = left + right

    # extrapolate
    # since we have repeated end-knots, only the last 2 basis functions are
    # non-zero at the end-knots, and they have equal and opposite gradient.
    if (any(x_extrapolte_r) or any(x_extrapolte_l)):
        bases[~x_interpolate] = 0.0

        denom = aug_knots[spline_order:-1] - aug_knots[: -spline_order - 1]
        left = prev_bases[:, :-1] / denom

        denom = aug_knots[spline_order + 1 :] - aug_knots[1:-spline_order]
        right = prev_bases[:, 1:] / denom

        grads = (spline_order) * (left - right)

        if any(x_extrapolte_l):
            val = grads[0] * x[x_extrapolte_l] + bases[-2]
            bases[x_extrapolte_l] = val
        if any(x_extrapolte_r):
            val = grads[1] * (x[x_extrapolte_r] - 1) + bases[-1]
            bases[x_extrapolte_r] = val
    # get rid of the added values at 0, and 1
    bases = bases[:-2]
    return bases

def calc_uncertainty(d,par):
    # Get the absolute path of the directory containing this file
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the path to LUT_par_uncertainty.txt
    lut_uncertainty_path = os.path.join(current_dir, 'LUT', 'LUT_par_uncertainty.txt')

    # Load the data
    f=open(lut_uncertainty_path,'r')
    lines=f.readlines()
    f.close()
    depth=[]
    bins=[]
    rms=[]
    
    num=int(len(lines)/3)
    for i in range(num):
        line1=lines[i*3]
        line2=lines[i*3+1]
        line3=lines[i*3+2]

        depth.append(int(line1))
        bins.append([float(i) for i in line2.split(',')])
        rms.append([float(i) for i in line3.split(',')])
    
    bins=np.array(bins)
    rms=np.array(rms)
    
    e=np.zeros(d.size)+np.nan

    
    for k in range(d.size):
        idx1=np.searchsorted(depth,d[k])
    
        outside=((idx1==0) | (idx1==71))
        toolow=(idx1==0)
    
        idx1=np.where(toolow,idx1,idx1-1)
        idx2=np.where(outside,idx1,idx1+1)

        i1=np.searchsorted(bins[idx1,:],par[k])
    
        toohigh= (i1==11)
        toolow= (i1==0)
    
        i1=np.where(toohigh,i1-1,i1)
    
        i2=np.where(toolow,i1,i1-1)
    
        e1=rms[idx1,i2]
    
    
        j1=np.searchsorted(bins[idx2,:],par[k])
    
        toohigh= (j1==11)
        toolow= (j1==0)
    
        j1=np.where(toohigh,j1-1,j1)
    
        j2=np.where(toolow,j1,j1-1)
    
        e2=rms[idx2,j2]     
    
        # slope=np.where(idx1==idx2,0,(e2-e1)/(depth[idx2]-depth[idx1])) #error Message when depth > 200 dbars
        # Modified from J. Tan’s version, based on her recommendations.
        if idx1 == idx2:
            slope = 0
        else:
            slope = (e2 - e1) / (depth[idx2] - depth[idx1])

        e[k]=slope*(d[k]-depth[idx1])+e1
    
    return e


def JTan_PAR_from_Ed(profile):
    f0_par=2411.2579039931697 
    f0_Ed=np.array([117.1379,195.4065,202.604,188.264])
    
    # profile is a 2-d array, n depth x 5 columns, i.e., [Ed380,Ed443,Ed490,Ed555,depth]
    # Get the absolute path of the directory containing this file
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the path to coeff.txt
    coeff_path = os.path.join(current_dir, 'LUT', 'coeff.txt')

    # Load the coefficients
    coef=np.genfromtxt(coeff_path,delimiter=',')
    coef=coef.reshape((100,4,100))

    b=b_spline_basis(profile[:,-1])
    
    par=np.zeros((profile.shape[0],100))+np.nan
    
    for j in range(100):
        par[:,j]=(np.dot(b,coef[j,0,:])*profile[:,0]/f0_Ed[0]+np.dot(b,coef[j,1,:])*profile[:,1]/f0_Ed[1]+\
        np.dot(b,coef[j,2,:])*profile[:,2]/f0_Ed[2]+np.dot(b,coef[j,3,:])*profile[:,3]/f0_Ed[3])*f0_par


    e1=np.std(par,1)

    e2=calc_uncertainty(profile[:,-1],np.mean(par,1))

    e=np.sqrt(e1**2+e2**2)

    return np.mean(par,1),e

# read example data, 1st column - depth, 2nd - PAR (micro E/m2/s)
# 3rd - 6th columns: Ed380 (mW/cm2/micron), Ed443 (mW/cm2/micron), Ed490 (mW/cm2/micron), Ed555 (mW/cm2/micron)

# df=read_csv('example/example_data.csv',skiprows=None)
# df=df.values
#
# d=df[:,0]
# par=df[:,1]
# profile=df[:,[2,3,4,5,0]]
#
# par_est,uncertainty=JTan_PAR_from_Ed(profile)
#
# print('Estimated PAR:')
# print(par_est)
#
# plt.figure(figsize=(10,4))
# plt.subplot(1,2,1)
# plt.plot(d[np.argsort(d)],par_est[np.argsort(d)],'b',label='modeled')
# plt.plot(d[np.argsort(d)],par[np.argsort(d)],'r',label='measured')
# plt.xlabel('Depth, m',fontsize=14)
# plt.ylabel('PAR, '+r'$\mu$E/m$^{2}$/s',fontsize=14)
# plt.tick_params(labelsize=12)
# plt.yscale('log')
# plt.legend(fontsize=12)
# plt.subplot(1,2,2)
# plt.plot(d[np.argsort(d)],uncertainty[np.argsort(d)],'k')
# plt.xlabel('Depth, m',fontsize=14)
# plt.ylabel(r'$\Delta$'+'PAR, '+r'$\mu$E/m$^{2}$/s',fontsize=14)
# plt.tick_params(labelsize=12)
# plt.tight_layout()
# plt.savefig('example/example_measured_vs_modeled.png')
# plt.close()
