import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from numba import jit
import scipy

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['ps.fonttype']=42
matplotlib.rcParams['axes.labelsize']=14
matplotlib.rcParams['xtick.labelsize']=12
matplotlib.rcParams['ytick.labelsize']=12
matplotlib.rcParams['font.sans-serif']="Arial"
matplotlib.rcParams['font.family']="sans-serif"

@jit
def calculate_fpt(kon, koff, klife, t):
    klt = klife
    psuc = koff*kon*(klt+kon)**(-1)*\
           (klt+2*koff+kon)*(klt**2+3*klt*koff+2*koff**2+2*klt*kon+koff*kon+kon**2)**(-1)

    fpt = (-1/2)*np.exp(1)**((-1/2)*(2*klt+3*koff+2*kon+koff**(1/2)*(koff+8*kon)**(1/2))*t)*koff**(1/2)*kon*(koff**2+7*koff*kon+(-8)*kon**2)**(-1)*((1+np.exp(1)**(koff**(1/2)*(koff+8*kon)**(1/2)*t)+(-2)*np.exp(1)**((1/2)*(3*koff*t+koff**(1/2)*(koff+8*kon)**(1/2)*t)))*koff**(3/2)+8*(1+np.exp(1)**(koff**(1/2)*(koff+8*kon)**(1/2)*t)+(-2)*np.exp(1)**((1/2)*(3*koff*t+koff**(1/2)*(koff+8*kon)**(1/2)*t)))*koff**(1/2)*kon+((-1)+np.exp(1)**(koff**(1/2)*(koff+8*kon)**(1/2)*t))*koff*(koff+8*kon)**(1/2)+2*((-1)+np.exp(1)**(koff**(1/2)*(koff+8*kon)**(1/2)*t))*kon*(koff+8*kon)**(1/2))/psuc

    return fpt

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

    # Get FPT, normalized so that the integral is 1.
    # I evaluate at the 'midpoints' of the empirical histogram
    t_FPT = np.arange(0.01, 3.98, 0.02)
    FPTout = calculate_fpt(kon, koff, klife, t_FPT)
    FPTout = FPTout / (np.sum(FPTout)*0.02)

    return tcout, csout, FPTout

@jit
def simulate_distributions_return_f1frac(kon, koff, klife, clus_size=10, nlim=2, nreals=200000):
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
    f1_frac = mm/nreals

    # Get FPT
    # I evaluate at the 'midpoints' of the empirical histogram
    t_FPT = np.arange(0.01, 3.98, 0.02)
    FPTout = calculate_fpt(kon, koff, klife, t_FPT)
    # Normalize FPT so AUC is 1
    FPTout = FPTout / (np.sum(FPTout)*0.02)

    return tcout, csout, FPTout, f1_frac


def visualize_distributions(kons, koffs, klifes, Ntagfile, tagLifeFile, FPTFile):

    # Load data
    # df_Ntag = pd.read_csv('N_tag.csv')
    # df_life = pd.read_csv('life_tag.csv')
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

    # Generate distributions
    Ntag_heights, _ = np.histogram(Ntag, bins=np.arange(10, 50), density=True)
    lifetag_heights, _ = np.histogram(lifetag, bins=np.arange(0, 50), density=True)
    FPTtag_heights, _ = np.histogram(FPTtag, bins=np.arange(0.01, 4., 0.02), density=True)


    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(8, 3.5))

    for kon, koff, klife in zip(kons, koffs, klifes):
        # Conduct simulation
        tc_sim, cs_sim, fpt_heights_sim = simulate_distributions(kon, koff, klife)
        Ntag_heights_sim, _ = np.histogram(cs_sim, bins=np.arange(10, 50), density=True)
        lifetag_heights_sim, _ = np.histogram(tc_sim, bins=np.arange(0, 50), density=True)

        ax[0].plot(np.arange(10.5, 49.5, 1), Ntag_heights_sim, label='sim', alpha=0.5, color='r')
        ax[1].plot(np.arange(0.5, 49.5, 1), lifetag_heights_sim, label='sim', alpha=0.5, color='r')
        ax[2].plot(np.arange(0, 3.98, 0.02), fpt_heights_sim, label='sim', alpha=0.5, color='r')



    ax[0].plot(np.arange(10.5, 49.5, 1), Ntag_heights, label='exp', color='k', linewidth=3)
    ax[1].plot(np.arange(0.5, 49.5, 1), lifetag_heights, label='exp', color='k', linewidth=3)
    ax[2].plot(np.arange(0, 3.98, 0.02), FPTtag_heights, label='sim', color='k', linewidth=3)

    ax[0].set_xlabel('tracks in cluster')
    ax[0].set_ylabel('probability density')

    ax[1].set_xlabel('cluster lifetime (s)')
    ax[1].set_ylabel('probability density')

    ax[2].set_xlabel('first passage time (s)')
    ax[2].set_ylabel('probability density')
    plt.show()

def find_nearest(array, value, range):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return range[np.argwhere(array == array[idx])[0][0]]

def get_credible_region(df, alpha):
    # Alpha = 0.05 = 95% credible, measure at 0.025 and 0.975 quartile
    # Alpha = 0.01 = 99% credible, measure at 0.005 and 0.995 quartile.
    df['weight'] = 1/df['error']
    kon_range = np.arange(0, 5, 0.01)
    koff_range = np.arange(0, 30, 0.01)
    klife_range = np.arange(0, 1, 0.01)

    kon_kde = scipy.stats.gaussian_kde(df['kon'].values,
                                       weights=df['weight'].values)
    koff_kde = scipy.stats.gaussian_kde(df['koff'].values,
                                        weights=df['weight'].values)
    klife_kde = scipy.stats.gaussian_kde(df['klife'].values,
                                         weights=df['weight'].values)

    kon_cdf = np.cumsum(kon_kde.evaluate(kon_range))
    kon_cdf = kon_cdf/np.max(kon_cdf)

    koff_cdf = np.cumsum(koff_kde.evaluate(koff_range))
    koff_cdf = koff_cdf/np.max(koff_cdf)

    klife_cdf = np.cumsum(klife_kde.evaluate(klife_range))
    klife_cdf = klife_cdf/np.max(klife_cdf)

    lq = alpha / 2
    uq = 1 - (alpha/2)

    # CDF values naturally represent quantiles [0,1]
    kon_lower = find_nearest(kon_cdf, lq, kon_range)
    kon_upper = find_nearest(kon_cdf, uq, kon_range)

    koff_lower = find_nearest(koff_cdf, lq, koff_range)
    koff_upper = find_nearest(koff_cdf, uq, koff_range)

    klife_lower = find_nearest(klife_cdf, lq, klife_range)
    klife_upper = find_nearest(klife_cdf, uq, klife_range)

    return kon_lower, kon_upper, koff_lower, koff_upper, klife_lower, klife_upper

def process_string_to_array(inputdf, key):
    nentries = len(inputdf)

    unpacked = []
    test = inputdf[key][1:-1]
    unpacked = np.array([x.lstrip() for x in test.split(',')], dtype=np.float64)


    return unpacked

def generate_report(df):
    df_best1 = df[((df['colocaliz_UB']+df['colocaliz_LB'])/2 <= 0.01)]
    df_best1 = df_best1[(df_best1['error'] <= np.percentile(df_best1['error'], 10))]

    df_best2 = df[((df['colocaliz_UB']+df['colocaliz_LB'])/2 <= 0.01)]
    df_best2 = df_best2[(df_best2['error'] <= np.percentile(df_best2['error'], 5))]

    df_best3 = df[((df['colocaliz_UB']+df['colocaliz_LB'])/2 <= 0.01)]
    df_best3 = df_best3[(df_best3['error'] <= np.percentile(df_best3['error'], 1))]

    # Compute the distribution of f1_frac for the three regions of parameter space
    df_best1_f1frac = []
    df_best2_f1frac = []
    df_best3_f1frac = []

    for _, entry in df_best1.iterrows():
        kon = entry['kon']
        koff = entry['koff']
        klife = entry['klife']

        _, _, _, f1frac = simulate_distributions_return_f1frac(kon, koff, klife)
        df_best1_f1frac.append(f1frac)
    for _, entry in df_best2.iterrows():
        kon = entry['kon']
        koff = entry['koff']
        klife = entry['klife']

        _, _, _, f1frac = simulate_distributions_return_f1frac(kon, koff, klife)
        df_best2_f1frac.append(f1frac)
    for _, entry in df_best3.iterrows():
        kon = entry['kon']
        koff = entry['koff']
        klife = entry['klife']

        _, _, _, f1frac = simulate_distributions_return_f1frac(kon, koff, klife)
        df_best3_f1frac.append(f1frac)

    df_best1_f1unfilt = df_best1.copy()
    df_best2_f1unfilt = df_best2.copy()
    df_best3_f1unfilt = df_best3.copy()

    df_best1 = df_best1[(0.4 < np.array(df_best1_f1frac)) & (np.array(df_best1_f1frac) < 0.8)]
    df_best2 = df_best2[(0.4 < np.array(df_best2_f1frac)) & (np.array(df_best2_f1frac) < 0.8)]
    df_best3 = df_best3[(0.4 < np.array(df_best3_f1frac)) & (np.array(df_best3_f1frac) < 0.8)]


    color1 = 'red'
    color2 = '#2176FF'
    color3 = '#A5C882'

    kon_range = [0, 3]
    koff_range = [0, 20]
    klife_range =[0, 0.3]

    kon_ticks = [0, 1, 2, 3]
    koff_ticks = [0, 10, 20]
    klife_ticks = [0, 0.1, 0.2, 0.3]

    kon_mean1 = np.mean(df_best1['kon'])
    koff_mean1 = np.mean(df_best1['koff'])
    klife_mean1 = np.mean(df_best1['klife'])
    kon_std1 = np.std(df_best1['kon'])
    koff_std1 = np.std(df_best1['koff'])
    klife_std1 = np.std(df_best1['klife'])

    kon_mean2 = np.mean(df_best2['kon'])
    koff_mean2 = np.mean(df_best2['koff'])
    klife_mean2 = np.mean(df_best2['klife'])
    kon_std2 = np.std(df_best2['kon'])
    koff_std2 = np.std(df_best2['koff'])
    klife_std2 = np.std(df_best2['klife'])

    kon_mean3 = np.mean(df_best3['kon'])
    koff_mean3 = np.mean(df_best3['koff'])
    klife_mean3 = np.mean(df_best3['klife'])
    kon_std3 = np.std(df_best3['kon'])
    koff_std3 = np.std(df_best3['koff'])
    klife_std3 = np.std(df_best3['klife'])

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(8, 6))
    ax[0, 0].plot(df_best1['kon'].values, df_best1['koff'].values, marker='s', color=color1, label='Top 10%', linestyle='None', markerfacecolor=color1)
    ax[0, 0].plot(df_best2['kon'].values, df_best2['koff'].values, marker='s', color=color2, label='Top 5%', linestyle='None', markerfacecolor=color2)
    ax[0, 0].plot(df_best3['kon'].values, df_best3['koff'].values, marker='s', color=color3, label='Top 1%', linestyle='None', markerfacecolor=color3)
    ax[0, 0].set_xlim(kon_range)
    ax[0, 0].set_ylim(koff_range)
    ax[0, 0].set_xticks(kon_ticks)
    ax[0, 0].set_yticks(koff_ticks)
    ax[0, 0].set_xlabel('$k_{on}$')
    ax[0, 0].set_ylabel('$k_{off}$')
    ax[0, 0].legend(loc='upper center', ncol=3, bbox_to_anchor=(1.2, 1.2))#, borderaxespad=0.)


    ax[0, 1].plot(df_best1['koff'].values, df_best1['klife'].values, marker='s', color=color1, linestyle='None', markerfacecolor=color1)
    ax[0, 1].plot(df_best2['koff'].values, df_best2['klife'].values, marker='s', color=color2, linestyle='None', markerfacecolor=color2)
    ax[0, 1].plot(df_best3['koff'].values, df_best3['klife'].values, marker='s', color=color3, linestyle='None', markerfacecolor=color3)
    ax[0, 1].set_xlim(koff_range)
    ax[0, 1].set_ylim(klife_range)
    ax[0, 1].set_xticks(koff_ticks)
    ax[0, 1].set_yticks(klife_ticks)
    ax[0, 1].set_xlabel('$k_{off}$')
    ax[0, 1].set_ylabel('$k_{lt}$')

    ax[0, 2].plot(df_best1['klife'].values, df_best1['kon'].values, marker='s', color=color1, linestyle='None', markerfacecolor=color1)
    ax[0, 2].plot(df_best2['klife'].values, df_best2['kon'].values, marker='s', color=color2, linestyle='None', markerfacecolor=color2)
    ax[0, 2].plot(df_best3['klife'].values, df_best3['kon'].values, marker='s', color=color3, linestyle='None', markerfacecolor=color3)
    ax[0, 2].set_xlim(klife_range)
    ax[0, 2].set_ylim(kon_range)
    ax[0, 2].set_xticks(klife_ticks)
    ax[0, 2].set_yticks(kon_ticks)
    ax[0, 2].set_xlabel('$k_{lt}$')
    ax[0, 2].set_ylabel('$k_{on}$')

    sns.kdeplot(data=df_best1, x='kon', ax=ax[1, 0], color=color1)
    sns.kdeplot(data=df_best2, x='kon', ax=ax[1, 0], color=color2)
    sns.kdeplot(data=df_best3, x='kon', ax=ax[1, 0], color=color3)
    ax[1, 0].set_xlabel('$k_{on}$')
    ax[1, 0].set_xlim(kon_range)
    # ax[1, 0].set_title('Top 1%: {:.1f}$\pm${:.2f}'.format(kon_mean3, kon_std3))

    sns.kdeplot(data=df_best1, x='koff', ax=ax[1, 1], color=color1)
    sns.kdeplot(data=df_best2, x='koff', ax=ax[1, 1], color=color2)
    sns.kdeplot(data=df_best3, x='koff', ax=ax[1, 1], color=color3)
    ax[1, 1].set_xlabel('$k_{off}$')
    ax[1, 1].set_xlim(koff_range)
    # ax[1, 1].set_title('Top 1%: {:.1f}$\pm${:.2f}'.format(koff_mean3, koff_std3))

    sns.kdeplot(data=df_best1, x='klife', ax=ax[1, 2], color=color1)
    sns.kdeplot(data=df_best2, x='klife', ax=ax[1, 2], color=color2)
    sns.kdeplot(data=df_best3, x='klife', ax=ax[1, 2], color=color3)
    ax[1, 2].set_xlabel('$k_{lt}$')
    ax[1, 2].set_xlim(klife_range)
    # ax[1, 2].set_title('Top 1%: {:.1f}$\pm${:.2f}'.format(klife_mean3, klife_std3))

    # plt.tight_layout()
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=.5, hspace=.4)

    plt.savefig('20210710_result_averagepSuc.png', dpi=300)
    plt.savefig('20210710_result_averagepSuc.pdf', dpi=300, transparent=True)

def generate_report_add_f1frac(df):
    df = pd.read_csv('pilot_sweep11_addfpt_weightFPT0-1_addpSuc_addHistDistrs.csv')
    # df_best1 = df[(df['colocaliz_UB'] <= 0.015)]
    # df_best1 = df_best1[(df_best1['error'] <= np.percentile(df_best1['error'], 10))]

    df_best2 = df[(df['colocaliz_UB'] <= 0.015)]
    df_best2 = df_best2[(df_best2['error'] <= np.percentile(df_best2['error'], 5))]

    # df_best3 = df[(df['colocaliz_UB'] <= 0.015)]
    # df_best3 = df_best3[(df_best3['error'] <= np.percentile(df_best3['error'], 1))]

    color1 = 'red'
    color2 = '#2176FF'
    color3 = '#A5C882'

    kon_range = [0, 3]
    koff_range = [0, 20]
    klife_range =[0, 0.3]

    kon_ticks = [0, 1, 2, 3]
    koff_ticks = [0, 10, 20]
    klife_ticks = [0, 0.1, 0.2, 0.3]

    # Compute the distribution of f1_frac for the three regions of parameter space
    # df_best1_f1frac = []
    df_best2_f1frac = []
    # df_best3_f1frac = []

    # for _, entry in df_best1.iterrows():
    #     kon = entry['kon']
    #     koff = entry['koff']
    #     klife = entry['klife']
    #
    #     _, _, _, f1frac = simulate_distributions_return_f1frac(kon, koff, klife)
    #     df_best1_f1frac.append(f1frac)
    for _, entry in df_best2.iterrows():
        kon = entry['kon']
        koff = entry['koff']
        klife = entry['klife']

        _, _, _, f1frac = simulate_distributions_return_f1frac(kon, koff, klife)
        df_best2_f1frac.append(f1frac)
    # for _, entry in df_best3.iterrows():
    #     kon = entry['kon']
    #     koff = entry['koff']
    #     klife = entry['klife']
    #
    #     _, _, _, f1frac = simulate_distributions_return_f1frac(kon, koff, klife)
    #     df_best3_f1frac.append(f1frac)

    # df_best1_f1unfilt = df_best1.copy()
    df_best2_f1unfilt = df_best2.copy()
    # df_best3_f1unfilt = df_best3.copy()

    # df_best1 = df_best1[(0.4 < np.array(df_best1_f1frac)) & (np.array(df_best1_f1frac) < 0.8)]
    df_best2 = df_best2[(0.4 < np.array(df_best2_f1frac)) & (np.array(df_best2_f1frac) < 0.8)]
    # df_best3 = df_best3[(0.4 < np.array(df_best3_f1frac)) & (np.array(df_best3_f1frac) < 0.8)]

    # kon_mean1 = np.mean(df_best1['kon'])
    # koff_mean1 = np.mean(df_best1['koff'])
    # klife_mean1 = np.mean(df_best1['klife'])
    # kon_std1 = np.std(df_best1['kon'])
    # koff_std1 = np.std(df_best1['koff'])
    # klife_std1 = np.std(df_best1['klife'])

    kon_mean2 = np.mean(df_best2['kon'])
    koff_mean2 = np.mean(df_best2['koff'])
    klife_mean2 = np.mean(df_best2['klife'])
    kon_std2 = np.std(df_best2['kon'])
    koff_std2 = np.std(df_best2['koff'])
    klife_std2 = np.std(df_best2['klife'])

    # kon_mean3 = np.mean(df_best3['kon'])
    # koff_mean3 = np.mean(df_best3['koff'])
    # klife_mean3 = np.mean(df_best3['klife'])
    # kon_std3 = np.std(df_best3['kon'])
    # koff_std3 = np.std(df_best3['koff'])
    # klife_std3 = np.std(df_best3['klife'])

    # fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(8, 8))

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(8.5*1.5, 3*1.5))
    # ax[0, 0].plot(df_best1['kon'].values, df_best1['koff'].values, marker='s', color=color1, label='Top 10%', linestyle='None', markerfacecolor=color1)
    # ax[0, 0].plot(df_best2['kon'].values, df_best2['koff'].values, marker='s', color=color2, label='Top 5%', linestyle='None', markerfacecolor=color2)
    # ax[0, 0].plot(df_best3['kon'].values, df_best3['koff'].values, marker='s', color=color3, label='Top 1%', linestyle='None', markerfacecolor=color3)
    # ax[0, 0].set_xlim(kon_range)
    # ax[0, 0].set_ylim(koff_range)
    # ax[0, 0].set_xticks(kon_ticks)
    # ax[0, 0].set_yticks(koff_ticks)
    # ax[0, 0].set_xlabel('$k_{on}$')
    # ax[0, 0].set_ylabel('$k_{off}$')
    # # ax[0, 0].legend(loc='upper center', ncol=3, bbox_to_anchor=(1.3, 1.3))#, borderaxespad=0.)
    #
    #
    # ax[0, 1].plot(df_best1['koff'].values, df_best1['klife'].values, marker='s', color=color1, linestyle='None', markerfacecolor=color1)
    # ax[0, 1].plot(df_best2['koff'].values, df_best2['klife'].values, marker='s', color=color2, linestyle='None', markerfacecolor=color2)
    # ax[0, 1].plot(df_best3['koff'].values, df_best3['klife'].values, marker='s', color=color3, linestyle='None', markerfacecolor=color3)
    # ax[0, 1].set_xlim(koff_range)
    # ax[0, 1].set_ylim(klife_range)
    # ax[0, 1].set_xticks(koff_ticks)
    # ax[0, 1].set_yticks(klife_ticks)
    # ax[0, 1].set_xlabel('$k_{off}$')
    # ax[0, 1].set_ylabel('$k_{life}$')
    #
    # ax[0, 2].plot(df_best1['klife'].values, df_best1['kon'].values, marker='s', color=color1, linestyle='None', markerfacecolor=color1)
    # ax[0, 2].plot(df_best2['klife'].values, df_best2['kon'].values, marker='s', color=color2, linestyle='None', markerfacecolor=color2)
    # ax[0, 2].plot(df_best3['klife'].values, df_best3['kon'].values, marker='s', color=color3, linestyle='None', markerfacecolor=color3)
    # ax[0, 2].set_xlim(klife_range)
    # ax[0, 2].set_ylim(kon_range)
    # ax[0, 2].set_xticks(klife_ticks)
    # ax[0, 2].set_yticks(kon_ticks)
    # ax[0, 2].set_xlabel('$k_{life}$')
    # ax[0, 2].set_ylabel('$k_{on}$')

    # sns.kdeplot(data=df_best1, x='kon', ax=ax[0], color=color1)
    sns.kdeplot(data=df_best2, x='kon', ax=ax[0], color=color2)
    # sns.kdeplot(data=df_best3, x='kon', ax=ax[0], color=color3)
    ax[0].set_xlabel('$k_{on}$')
    ax[0].set_xlim(kon_range)
    ax[0].set_title('Top 5%: {:.1f}$\pm${:.2f}'.format(kon_mean2, kon_std2))

    # sns.kdeplot(data=df_best1, x='koff', ax=ax[1], color=color1)
    sns.kdeplot(data=df_best2, x='koff', ax=ax[1], color=color2)
    # sns.kdeplot(data=df_best3, x='koff', ax=ax[1], color=color3)
    ax[1].set_xlabel('$k_{off}$')
    ax[1].set_xlim(koff_range)
    ax[1].set_title('Top 5%: {:.1f}$\pm${:.2f}'.format(koff_mean2, koff_std2))

    # sns.kdeplot(data=df_best1, x='klife', ax=ax[2], color=color1)
    sns.kdeplot(data=df_best2, x='klife', ax=ax[2], color=color2)
    # sns.kdeplot(data=df_best3, x='klife', ax=ax[2], color=color3)
    ax[2].set_xlabel('$k_{life}$')
    ax[2].set_xlim(klife_range)
    ax[2].set_title('Top 5%: {:.2f}$\pm${:.3f}'.format(klife_mean2, klife_std2))


    # ax[2, 0].plot(df_best1_f1unfilt['kon'].values, df_best1_f1frac, marker='s', color=color1, linestyle='None', markerfacecolor=color1)
    # ax[2, 0].plot(df_best2_f1unfilt['kon'].values, df_best2_f1frac, marker='s', color=color2, linestyle='None', markerfacecolor=color2)
    # ax[2, 0].plot(df_best3_f1unfilt['kon'].values, df_best3_f1frac, marker='s', color=color3, linestyle='None', markerfacecolor=color3)
    # ax[2, 0].set_xlim(kon_range)
    # ax[2, 0].set_xticks(kon_ticks)
    # ax[2, 0].set_xlabel('$k_{on}$')
    # ax[2, 0].set_ylabel('$f_1$')
    #
    # ax[2, 1].plot(df_best1_f1unfilt['koff'].values, df_best1_f1frac, marker='s', color=color1, linestyle='None',
    #               markerfacecolor=color1)
    # ax[2, 1].plot(df_best2_f1unfilt['koff'].values, df_best2_f1frac, marker='s', color=color2, linestyle='None',
    #               markerfacecolor=color2)
    # ax[2, 1].plot(df_best3_f1unfilt['koff'].values, df_best3_f1frac, marker='s', color=color3, linestyle='None',
    #               markerfacecolor=color3)
    # ax[2, 1].set_xlim(koff_range)
    # ax[2, 1].set_xticks(koff_ticks)
    # ax[2, 1].set_xlabel('$k_{off}$')
    # ax[2, 1].set_ylabel('$f_1$')
    #
    # ax[2, 2].plot(df_best1_f1unfilt['klife'].values, df_best1_f1frac, marker='s', color=color1, linestyle='None',
    #               markerfacecolor=color1)
    # ax[2, 2].plot(df_best2_f1unfilt['klife'].values, df_best2_f1frac, marker='s', color=color2, linestyle='None',
    #               markerfacecolor=color2)
    # ax[2, 2].plot(df_best3_f1unfilt['klife'].values, df_best3_f1frac, marker='s', color=color3, linestyle='None',
    #               markerfacecolor=color3)
    # ax[2, 2].set_xlim(klife_range)
    # ax[2, 2].set_xticks(klife_ticks)
    # ax[2, 2].set_xlabel('$k_{life}$')
    # ax[2, 2].set_ylabel('$f_1$')


    # print('test')
    plt.tight_layout()
    # plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=.5, hspace=.5)
    # plt.show()

    # ax[0, 0].plot(2.2, 5.1, 'kx')
    # ax[0, 1].plot(5.1, 0.049, 'kx')
    # ax[0, 2].plot(0.049, 2.2, 'kx')

    plt.savefig('20210702_result_v3_pFail0-15_show2-2_5-1_0-049.png', dpi=300)
    plt.savefig('20210702_result_v3_pFail0-15_show2-2_5-1_0-049.pdf', dpi=300, transparent=True)


def __main__():
    df = pd.read_csv('pilot_sweep11_addfpt_weightFPT0-1_addpSuc_addHistDistrs.csv')
    # generate_report_add_f1frac(df)
    generate_report(df)

if __name__ == '__main__':
    __main__()
