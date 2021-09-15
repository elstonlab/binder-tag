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
                          'initial_review')

def load_df():
    datapath = os.path.join(local_path,'sizes_nMinTrack.csv')
    df=pd.read_csv(datapath)
    df=np.sqrt(df/np.pi)
    df.columns=['tS-4', 'B-4',
                'tS-5', 'B-5',
                'tS-6', 'B-6',
                'tS-7', 'B-7',
                'tS-8', 'B-8',
                'tS-9', 'B-9',
                'tS-10', 'B-10',
                'tS-11', 'B-11',
                'tS-12', 'B-12',
                'tS-13', 'B-13',
                'tS-14', 'B-14',
                'tS-15', 'B-15',
                'tS-16', 'B-16',
                'tS-17', 'B-17',
                'tS-18', 'B-18',
                'tS-19', 'B-19',
                'tS-20', 'B-20']
    df=df.melt()
    return df

df=load_df()
plt.figure(figsize=(5.5,4))
sns.violinplot(x='variable',y='value',data=df,bw=1,cut=0,inner='box')
plt.ylabel('radius (nm)')
plt.xlabel('Minimum tracks to define clusters (tagSrc or Binder)')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(local_path,'nMinTrack_clusterRadius_small_20210902.png'),dpi=300,bbox_inches = "tight")
plt.savefig(os.path.join(local_path,'nMinTrack_clusterRadius_small_20210902.pdf'),dpi=300,bbox_inches = "tight",transparent=True)

#plt.xlim(50,125)
#plt.xticks(np.arange(50,150,25))
#plt.tight_layout()
#plt.savefig(os.path.join(local_path,'bootstrapped_radii.png'),dpi=300,bbox_inches = "tight")
#plt.savefig(os.path.join(local_path,'bootstrapped_radii.pdf'),dpi=300,bbox_inches = "tight",transparent=True)
