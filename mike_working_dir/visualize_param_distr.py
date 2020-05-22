import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np

import matplotlib
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['ps.fonttype']=42
matplotlib.rcParams['font.sans-serif']="Arial"
matplotlib.rcParams['font.family']="sans-serif"

evolved_params=pd.read_csv('mike_working_dir\\UpdatedIndividuals.csv')

evolved_rates = evolved_params[['k1','k2','k3','k4','k6','k7','k8']]
evolved_f1 = evolved_params[['f1']]
evolved_errors=evolved_params[['error']]

# The parameters were allowed run range bewteen 10-4 and 20s-1, except for "f1",
# which was allowed to range between 10-4 and 0.76.

# There are 28 evolved individuals, which came from 28 independent evolutions.
# Each evolution contained 100 individuals and ran for 40 generations, with
# a mutation rate of 0.1 and crossover rate of 0.5

fig=plt.figure()
gridspec.GridSpec(1,4)

plt.subplot2grid((1,4),(0,0),colspan=3,rowspan=1)
sns.swarmplot(data=evolved_rates.melt(),x='variable',y='value')
plt.ylabel('Parameter value')
plt.plot([0,1,2,3,4,5,6],np.array(evolved_rates.iloc[0]),'k',marker='D',markerfacecolor='white',linestyle='None',zorder=100)
plt.ylim(-.5,20.5)
plt.xlabel('')

plt.subplot2grid((1,4),(0,3),colspan=1,rowspan=1)
sns.swarmplot(data=evolved_f1.melt(),x='variable',y='value')
plt.plot([0],np.array(evolved_f1.iloc[0]),'k',marker='D',markerfacecolor='white',linestyle='None',zorder=100)
plt.ylabel('Parameter value')
plt.ylim(-.05,1.05)
plt.xlabel('')

plt.tight_layout()

plt.savefig('mike_working_dir\\evolved_parameters.png',dpi=300,bbox_inches='tight')
plt.savefig('mike_working_dir\\evolved_parameters.pdf',dpi=300,bbox_inches='tight')

plt.figure()
sns.swarmplot(data=evolved_errors.melt(),x='variable',y='value')
plt.plot([0],np.array(evolved_errors.iloc[0]),'k',marker='D',markerfacecolor='white',linestyle='None',zorder=100)
plt.ylabel('Error')
plt.xlabel('')

plt.savefig('mike_working_dir\\errors.png',dpi=300,bbox_inches='tight')
plt.savefig('mike_working_dir\\errors.pdf',dpi=300,bbox_inches='tight')
