import pandas as pd
import matplotlib.pyplot as plt
test=pd.read_csv('mike_working_dir/UpdatedIndividuals.csv')

plt.hist(test['k1'])
plt.hist(test['k2'])
plt.hist(test['k3'])
plt.hist(test['k4'])
plt.hist(test['k6'])
plt.hist(test['k7'])
plt.hist(test['k8'])
