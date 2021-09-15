import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['ps.fonttype']=42
matplotlib.rcParams['font.sans-serif']="Arial"
matplotlib.rcParams['font.family']="sans-serif"
local_path = os.path.join('C:\\','Users','mike','Documents','GitHub','binder-tag',
                          'initial_review','CAtagSrc_vs_Binder_analysis',
                          '20201212','followupanalysis','sizes_slow',
                          'bootstrapped')

def load_df():
    datapath = os.path.join(local_path,'slowTS-slowB-fastTS-fastB_ar50_bootstrapped.csv')
    df=pd.read_csv(datapath)
    df=np.sqrt(df/np.pi)
    df.columns=['tagSrc-CA, slow','Binder, slow','tagSrc-CA, fast','Binder, fast']
    df=df.melt()
    return df

df=load_df()
plt.figure(figsize=(2.5,2))
sns.violinplot(x='value',y='variable',data=df,bw=1,cut=0,inner='box')
plt.xlabel('radius (nm)')
plt.ylabel('')
# plt.xlim(75,150)
# plt.xticks(np.arange(75,175,25))
plt.xlim(50,125)
plt.xticks(np.arange(50,150,25))
plt.tight_layout()
plt.savefig(os.path.join(local_path,'bootstrapped_radii.png'),dpi=300,bbox_inches = "tight")
plt.savefig(os.path.join(local_path,'bootstrapped_radii.pdf'),dpi=300,bbox_inches = "tight",transparent=True)
