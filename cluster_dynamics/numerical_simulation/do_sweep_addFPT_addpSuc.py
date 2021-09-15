import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from numba import jit

from pymcmcstat import MCMC

@jit
def calculate_fpt(kon, koff, klife, t):
    klt = klife
    psuc = koff*kon*(klt+kon)**(-1)*\
           (klt+2*koff+kon)*(klt**2+3*klt*koff+2*koff**2+2*klt*kon+koff*kon+kon**2)**(-1)

    fpt = (-1/2)*np.exp(1)**((-1/2)*(2*klt+3*koff+2*kon+koff**(1/2)*(koff+8*kon)**(1/2))*t)*koff**(1/2)*kon*(koff**2+7*koff*kon+(-8)*kon**2)**(-1)*((1+np.exp(1)**(koff**(1/2)*(koff+8*kon)**(1/2)*t)+(-2)*np.exp(1)**((1/2)*(3*koff*t+koff**(1/2)*(koff+8*kon)**(1/2)*t)))*koff**(3/2)+8*(1+np.exp(1)**(koff**(1/2)*(koff+8*kon)**(1/2)*t)+(-2)*np.exp(1)**((1/2)*(3*koff*t+koff**(1/2)*(koff+8*kon)**(1/2)*t)))*koff**(1/2)*kon+((-1)+np.exp(1)**(koff**(1/2)*(koff+8*kon)**(1/2)*t))*koff*(koff+8*kon)**(1/2)+2*((-1)+np.exp(1)**(koff**(1/2)*(koff+8*kon)**(1/2)*t))*kon*(koff+8*kon)**(1/2))/psuc

    return fpt

@jit
def calculate_cluster_colocalization_probability_range(kon, koff, klife):
    # 5000, assumed number of observed tagsr clusters.
    # the effect of this metric on the number of clusters detected should
    # be robust to this constant
    NTag = 5000

    klt = klife
    psuc = klt*koff**10*kon**10*(klt+kon)**(-10)*(klt+2*koff+kon)**10* \
           (klt**2+3*klt*koff+2*koff**2+2*klt*kon+koff*kon+kon**2)**(-10)* \
           (klt*(klt+koff)*(klt+2*koff)+3*klt*(klt+koff)*kon+ \
            3*klt*kon**2+kon**3)**(-1)*(klt**2+2*koff**2+3*koff*kon+3*kon**2+3*klt*(koff+kon))
    # total number of clusters
    NTot = NTag/psuc

    #Binomial distribution for coincident clusters with probalitity of success = psuc^2
    NBT = NTot * psuc**2

    std_BT = NTot**(1/2) * (psuc**2 - psuc**4)**(1/2)



    # Our expectation, from our data:
    # (NBT+std_BT)/NTot < 0.01 (fewer than 1% of clusters with colocalization)
    lower = (NBT-std_BT)/NTot
    upper = (NBT+std_BT)/NTot
    return lower, upper

@jit
def simulate_distributions(kon, koff, klife, clus_size=10, nlim=2, nreals=200000):
    kk = 0
    mm = 0
    # tavg = 0

    # Preallocate, then drop nan values.
    csout = np.nan*np.ones(nreals)
    tcout = np.nan*np.ones(nreals)

    for i in range(nreals):
        n = 0
        m = 0
        flag = 0
        t = 0
        tc = 0
        while not flag:
            dt = -1 / (kon+koff*n+klife)*np.log(np.random.random())
            t = t+dt
            p = np.random.random()

            if kon/(kon+koff*n+klife) > p:
                n = n + 1
                if n > nlim:
                    flag = 1
                    mm = mm + 1
            elif ((kon+koff*n) / (kon+koff*n+klife)) > p:
                n = n - 1
                if n == 0:
                    m = m + 1
                    tc = t
            else:
                flag = 2
        if m >= clus_size and flag == 2:
            csout[i] = m
            tcout[i] = tc
            kk = kk + 1
            # tavg = (kk-1)*tavg/kk + tc/kk

    csout = csout[~np.isnan(csout)]
    tcout = tcout[~np.isnan(tcout)]
    # success_frac = kk/nreals
    # f1_frac = mm/nreals

    # Get FPT
    # I evaluate at the 'midpoints' of the empirical histogram
    t_FPT = np.arange(0.01, 3.98, 0.02)
    FPTout = calculate_fpt(kon, koff, klife, t_FPT)
    # Normalize FPT so AUC is 1
    FPTout = FPTout / (np.sum(FPTout)*0.02)

    return tcout, csout, FPTout


def do_score(kon, koff, klife, Ntag_heights, lifetag_heights, FPTtag_heights):
    # Conduct simulation
    try:
        tc_sim, cs_sim, FPT_heights_sim = simulate_distributions(kon, koff, klife)
        # # Normalize to maximum 1
        # FPT_heights_sim = FPT_heights_sim / np.max(FPT_heights_sim)

        if (not len(tc_sim)==0) and (not len(cs_sim)==0):
            # Generate distributions
            Ntag_heights_sim, c = np.histogram(cs_sim, bins=np.arange(0, 50), density=True)
            lifetag_heights_sim, d = np.histogram(tc_sim, bins=np.arange(0, 50), density=True)

            MSE_Ntag = np.mean(np.power(Ntag_heights[10:] - Ntag_heights_sim[10:], 2))
            MSE_lifetag = np.mean(np.power(lifetag_heights - lifetag_heights_sim, 2))

            # Empirical point for ignoring the first passage times.
            # Very short FPT are probably artifacts, since the cluster must exist
            # long enough for 10 tracks to arrive and leave.
            FPT_tmin = 20  # 20*0.02s = 0.4s
            MSE_FPTtag = np.mean(np.power(FPTtag_heights[FPT_tmin:] - FPT_heights_sim[FPT_tmin:], 2))

            # Empirically weight FPT at 1:10
            MSE_FPTtag = MSE_FPTtag*0.1
        else:
            MSE_Ntag = np.inf
            MSE_lifetag = np.inf
            MSE_FPTtag = np.inf
            Ntag_heights_sim = None
            lifetag_heights_sim = None
            FPT_heights_sim = None
    except:
        MSE_Ntag = np.inf
        MSE_lifetag = np.inf
        MSE_FPTtag = np.inf
        Ntag_heights_sim = None
        lifetag_heights_sim = None
        FPT_heights_sim = None

    return MSE_Ntag + MSE_lifetag + MSE_FPTtag, Ntag_heights_sim, lifetag_heights_sim, FPT_heights_sim

def do_sweep(Ntagfile, tagLifeFile, FPTFile, kons, koffs, klifes):
    # Load measured data
    df_Ntag = pd.read_csv(Ntagfile)
    df_life = pd.read_csv(tagLifeFile)
    df_FPT = pd.read_csv(FPTFile)

    Ntag = np.ravel(df_Ntag.values)
    lifetag = np.ravel(df_life.values)
    FPTtag = np.ravel(df_FPT.values)

    # Cleanup data (drop inf, nan)
    Ntag_tokeep = ~np.array(np.isnan(Ntag) | np.isinf(Ntag))
    lifetag_tokeep = ~np.array(np.isnan(lifetag) | np.isinf(lifetag))
    FPTtag_tokeep = ~np.array(np.isnan(FPTtag) | np.isinf(FPTtag))

    Ntag = Ntag[Ntag_tokeep]
    lifetag = lifetag[lifetag_tokeep]
    FPTtag = FPTtag[FPTtag_tokeep]


    # Precalculate the histograms for scoring. Note that we omit [0, 10) for scoring the Ntag_heights histograms
    Ntag_heights, _ = np.histogram(Ntag, bins=np.arange(0, 50), density=True)
    lifetag_heights, _ = np.histogram(lifetag, bins=np.arange(0, 50), density=True)
    FPTtag_heights, _ = np.histogram(FPTtag, bins=np.arange(0.01, 4, 0.02), density=True)

    # # Normalize to maximum 1
    df = pd.DataFrame(columns=['kon', 'koff', 'klife', 'error'])

    counter = 0
    total_parameter_sets = len(kons)*len(koffs)*len(klifes)
    for kon in kons:
        for koff in koffs:
            for klife in klifes:
                print('Scoring parameter set {:d} out of {:d}'.format(counter, total_parameter_sets))
                error, Ntag_heights_sim, lifetag_heights_sim, FPT_heights_sim = do_score(kon, koff, klife, Ntag_heights, lifetag_heights, FPTtag_heights)
                coloc_LB, coloc_UB = calculate_cluster_colocalization_probability_range(kon, koff, klife)
                df = df.append({'kon': kon, 'koff': koff, 'klife': klife,
                                'error': error, 'colocaliz_LB': coloc_LB,
                                'colocaliz_UB': coloc_UB,
                                'Ntag_heights_sim': Ntag_heights_sim.tolist(),
                                'lifetag_heights_sim': lifetag_heights_sim.tolist(),
                                'FPT_heights_sim': FPT_heights_sim.tolist()}, ignore_index=True)
                counter += 1

    return df


def __main__():
    kons = np.arange(0.9, 3.4, 0.05)
    koffs = np.arange(3, 21, 0.5)
    klifes = np.arange(0.01, 0.51, 0.01)


    df = do_sweep('N_tag.csv', 'life_tag.csv', 'all_Ti_tag.csv', kons, koffs, klifes)

    df.to_csv('pilot_sweep11_addfpt_weightFPT0-1_addpSuc_addHistDistrs.csv')

    errors = df['error'].values
    errors = np.reshape(errors, (len(kons),
                                 len(koffs),
                                 len(klifes)))
    errors = np.vstack((
                       np.hstack((errors[:, 0, :],
                                  errors[:, 1, :],
                                  errors[:, 2, :],
                                  errors[:, 3, :],
                                  errors[:, 4, :],
                                  errors[:, 5, :])),
                       np.hstack((errors[:, 6, :],
                                  errors[:, 7, :],
                                  errors[:, 8, :],
                                  errors[:, 9, :],
                                  errors[:, 10, :],
                                  errors[:, 11, :])),
                       np.hstack((errors[:, 12, :],
                                  errors[:, 13, :],
                                  errors[:, 14, :],
                                  errors[:, 15, :],
                                  errors[:, 16, :],
                                  errors[:, 17, :])),
                       np.hstack((errors[:, 18, :],
                                  errors[:, 19, :],
                                  errors[:, 20, :],
                                  errors[:, 21, :],
                                  errors[:, 22, :],
                                  errors[:, 23, :])),
                       np.hstack((errors[:, 24, :],
                                  errors[:, 25, :],
                                  errors[:, 26, :],
                                  errors[:, 27, :],
                                  errors[:, 28, :],
                                  errors[:, 29, :])),
                       np.hstack((errors[:, 30, :],
                                  errors[:, 31, :],
                                  errors[:, 32, :],
                                  errors[:, 33, :],
                                  errors[:, 34, :],
                                  errors[:, 35, :]))
    ))

    plt.figure()
    ax=plt.gca()
    sns.heatmap(np.log10(errors/np.nanmin(errors)), ax=ax)
    plt.show()

if __name__ == '__main__':
    __main__()
