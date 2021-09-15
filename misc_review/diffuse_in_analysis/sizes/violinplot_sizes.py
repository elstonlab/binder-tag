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

# C:\Users\mike\Documents\GitHub\binder-tag\initial_review\CAtagSrc_vs_Binder_analysis\20201212\followupanalysis\sizes\bootstrap_cluster_sizes
local_path = os.path.join('C:\\','Users','mike','Documents','GitHub','binder-tag',
                          'initial_review','diffuse_in_analysis',
                          'sizes',
                          'bootstrap_cluster_sizes')

def load_df():
    datapath = os.path.join(local_path,'tagsrc_binder_ar50_removeDiffuseIn_bootstrapped.csv')
    df=pd.read_csv(datapath)
    df=np.sqrt(df/np.pi)
    df.columns=['tagSrc','Binder']
    df=df.melt()
    return df

df=load_df()
plt.figure(figsize=(2.5,2))
sns.violinplot(x='value',y='variable',data=df,bw=1,cut=0,inner='box')
plt.xlabel('radius (nm)')
plt.title('Removed tracks diffusing in')
plt.ylabel('')
plt.xlim(25,125)
plt.xticks(np.arange(25,125,25))
plt.tight_layout()
plt.savefig(os.path.join(local_path,'bootstrapped_radii_removeDiffuseIn.png'),dpi=300,bbox_inches = "tight")
plt.savefig(os.path.join(local_path,'bootstrapped_radii_removeDiffuseIn.pdf'),dpi=300,bbox_inches = "tight",transparent=True)
