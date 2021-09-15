# import required packages
import numpy as np
import scipy.io as sio
import pandas as pd
import pickle
from datetime import datetime
from scipy.integrate import odeint
from pymcmcstat import MCMC, structures, plotting, propagation
import matlab.engine 
import matplotlib.pyplot as plt

np.seterr(over='ignore')

bt_data = sio.loadmat('BT_all_exp_data.mat',struct_as_record=True)

data_names = ['tag_codiff','bin_codiff','tag_all','cat','proc_1', 'proc_2', 'proc_3','proc_4']

exp_lft_data = []
exp_lft_x = bt_data['all_exp_data']['tag_codiff'][0][0]['X'][0][0][0]
exp_cat_data = []
exp_cat_x = bt_data['all_exp_data']['cat'][0][0]['X'][0][0][0]
exp_proc_data = []
exp_proc_x = bt_data['all_exp_data']['proc_1'][0][0]['X'][0][0][0]

exp_data = []
exp_x = []
for i,d in enumerate(data_names):
    exp_x.append(bt_data['all_exp_data'][d][0][0]['X'][0][0][0])
    exp_data.append(bt_data['all_exp_data'][d][0][0]['Y'][0][0][0])

#pad data with nan to be the same length/shape
pad = np.max([len(k) for k in exp_data])
exp_data = np.array([np.concatenate([k.flatten(), np.full(pad-len(k), np.nan)]) for k in exp_data])
exp_x =  np.array([np.concatenate([k.flatten(), np.full(pad-len(k), np.nan)]) for k in exp_x])
exp_data = np.array(exp_data).T
exp_x = np.array(exp_x.T)

eng = matlab.engine.start_matlab() #start matlab engine 

#get MSE between exp. and sim.
def BT_mse(params,mc_data):
    
    def mse(A,B):
        A = A[~np.isnan(A)]
        B = B[~np.isnan(B)]
        return(np.mean(np.subtract(A,B)**2))

    ndp, nbatch = mc_data.shape[0]
    exp_lft = np.array(mc_data.ydata).T[:3] 
    exp_cat = np.array(mc_data.ydata).T[3] 
    exp_proc = np.array(mc_data.ydata).T[4:] 

    sim_dists = run_BT(params)
    if np.shape(sim_dists) != (75,8): 
        return(np.array(sim_dists).T)
    sim_lft = np.array(sim_dists.T[:3])
    sim_cat =  np.array(sim_dists.T[3])
    sim_proc =  np.array(sim_dists.T[4:])

    mse_lft = [mse(exp_lft[i],sim_lft[i]) for i in range(3)]
    mse_cat = [mse(exp_cat,sim_cat)]
    mse_proc = [mse(exp_proc[i],sim_proc[i]) for i in range(4)]
    return(np.hstack([mse_lft,mse_cat,mse_proc]))

def run_BT(params, data = None):
    params = [float(p) for p in params]
    out = eng.markov_distributions(*params)
    out = [np.asarray(out[i]) for i in range(len(out))]
    pad = len(max(out, key=len))

    out = np.array([np.concatenate([d.flatten(), np.full(pad-len(d), np.nan)]) for d in out])
    
    return(np.array(out).T)


# initialize MCMC object
mcstat = MCMC.MCMC()

mcstat.data.add_data_set(x=np.arange(0,75),
                         y=exp_data,
                         user_defined_object=exp_x)


# add model parameters, initialize at best EA parameter set
mcstat.parameters.add_model_parameter(name='k1', theta0=5.31, minimum=0)
mcstat.parameters.add_model_parameter(name='k2', theta0=8.88, minimum=0)
mcstat.parameters.add_model_parameter(name='k3', theta0=1.39, minimum=0)
mcstat.parameters.add_model_parameter(name='k4', theta0=2.98, minimum=0)
mcstat.parameters.add_model_parameter(name='k6', theta0=0.63, minimum=0)
mcstat.parameters.add_model_parameter(name='k7', theta0=12.1, minimum=0)
mcstat.parameters.add_model_parameter(name='k8', theta0=0.23, minimum=0)
mcstat.parameters.add_model_parameter(name='f1', theta0=0.15, minimum=0,maximum=0.76)


# Generate options
mcstat.simulation_options.define_simulation_options(
    nsimu=1.0e2, updatesigma=True, 
    verbosity=False,save_to_json=True,
    save_lightly=False, waitbar=False )
     #save_to_json=True, verbosity=0, waitbar=True, save_to_bin=True)
# Define model object:
mcstat.model_settings.define_model_settings(
    sos_function=BT_mse,
    nbatch = 8,
    sigma2=0.01**2,S20=0.01*np.ones(8),N0=5*np.ones(8))

# Run simulation
mcstat.run_simulation()
# # Rerun starting from results of previous run
mcstat.simulation_options.nsimu = int(1.0e3)
mcstat.run_simulation(use_previous_results=True)

results = mcstat.simulation_results.results
chain = results['chain']
s2chain = results['s2chain']
names = results['names']

intervals = propagation.calculate_intervals(chain[500:], results, mcstat.data, run_BT,
                                                waitbar=False, nsample=250,s2chain = s2chain[500:])

now = datetime.now()
dt_string = now.strftime("%Y%m%d_%s")
save_file = dt_string + '_Intervals_S001N5.pickled'

pickle.dump(intervals, open(save_file, "wb" ) )

eng.quit()
