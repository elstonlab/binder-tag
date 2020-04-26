import os, errno
import numpy as np
import pickle
from operator import itemgetter
figsavedir = 'movie/'

try:
    os.makedirs(figsavedir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

def load_data(dir_to_check, score_cutoff):
    val_to_keep_end_score = score_cutoff

    #load new data:
    arr_best_scores = []
    arr_best_inds = []
    arr_end_scores = []
    files = os.listdir(dir_to_check)
    for i in range(0,len(files)):
        filename = dir_to_check + '/' + files[i]
        if os.path.isfile(filename):
            if '.pickled' in files[i]:
                arr_to_unpickle = pickle.load(open(filename,'rb'))
                arr_best_score, arr_best_ind = arr_to_unpickle
                temp_end_score = arr_best_score[-1]
                if temp_end_score < val_to_keep_end_score:
                    arr_end_scores.append(temp_end_score)
                    arr_best_scores.append(arr_best_score)
                    arr_best_inds.append(arr_best_ind)
        else:
            break

    arr_best_scores = np.asarray(arr_best_scores)
    arr_best_inds = np.asarray(arr_best_inds)
    arr_end_scores = np.asarray(arr_end_scores)

    print('Loaded ' + str(len(arr_best_scores)) + ' files (out of ' + str(len(files)) + ') with cutoff score of ' + str(score_cutoff))

    return arr_best_scores, arr_end_scores, arr_best_inds

################################################################
# Loads the data.
################################################################

arr_best_scores, arr_end_scores, arr_best_inds =\
load_data('mike_working_dir/BT_Updated_1001/',10000000)


np.set_printoptions(precision=1, suppress=True)

data = sorted(zip(arr_best_scores, arr_end_scores, arr_best_inds), key=itemgetter(1))

#MSE PLOT
import matplotlib.pyplot as plt
%matplotlib inline

gen = 40

range(1,targ_gen)

fig = plt.figure(figsize=(10,5))
for targ_gen in range(1,gen+1):
    for i in range(len(data)):
        plt.xlabel('Generation', fontsize=20, fontname='Arial')
        plt.ylabel('Mean square error', fontsize=20, fontname='Arial')
        plt.xlim([0,gen])
        plt.ylim([0,0.012])
        plt.plot(range(1,targ_gen+1), data[i][0][:targ_gen])
        plt.tick_params(axis='both',which='major',labelsize=14)
        cax = plt.gca()
        for tick in cax.get_xticklabels():
            tick.set_fontname('Arial')
        for tick in cax.get_yticklabels():
            tick.set_fontname('Arial')
    plt.savefig(figsavedir+'/MSE_timecourse_gen{:02d}.png'.format(targ_gen),dpi=100,bbox_inches='tight')
    plt.clf()

np.shape(data[0][0])
plt.plot(range(1,41),data[0][0][:(targ_gen+1)])
