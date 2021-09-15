#Binder-Tag Markov Chain Kinetics Model

# Markov chain model

##Source data:
- Lifetime distributions for Binder tracks in co-diffusion events are at `markov_chain_model/binder_lifetimes.csv`.
- Lifetime distributions for all Binder tracks (co-diffusion or not) are at `markov_chain_model/All events histogram.csv`
- Lifetime distributions for tagSrc tracks in co-diffusion events are at `markov_chain_model/tagsrc_lifetimes.csv`.
- Lifetime distributions for all tagSrc tracks (co-diffusion or not) are at `markov_chain_model/All events histogram.csv`
- Process times are at `markov_chain_model/event_times.mat`. This file has the same contents, though relabeled for convenience, as `show_dox_indept/checked_codiff_result_100pg.mat`.

##Simulation:
	The model had four possible states: state 1) closed tagSrc, state 2) open tagSrc, state 3 ) open tagSrc with Binder, and state 4) tagSrc dissociated from the plasma membrane. The Markov Chain model was simulated using a custom version of the dtmc() and simulate() functions. Simulations could begin in states 1, 2, or 3, and ended when state 4 was entered. For each starting state, 5000 Discrete-Time Markov-Chains were simulated. Similar to the experimental data, simulations were thrown out where state 3 was entered more than once. For simulations where state 3 was never entered, only the tagSrc lifetime was measured.

##Co-diffusion categories:
	There were four possible co-diffusion categories: category 1) tagSrc arrived first, Binder arrived then left, and finally tagSrc left, category 2) tagSrc arrived first, Binder arrived, and tagSrc and Binder left together, category 3) tagSrc and Binder arrived together, Binder left, then tagSrc left, and category 4) tagSrc and Binder both arrived and left together. Experimentally, in 76% of co-diffusion cases tagSrc arrived first (categories 1 and 2), and in 24% of cases tagSrc and Binder arrived together (categories 3 and 4). Thus, for the simulated co-diffusion categories, the proportion between categories 1 and 2 was measured and weighted by 0.76 and the proportion between categories 3 and 4 was measured and weighted by 0.24.

##Process distributions.
	There were four measurable processes we examined: process 1) time until tagSrc activation, process 2) time until tagSrc inactivation, process 3) time until Src dissociated from the plasma membrane and process 4) time until active Src dissociated from the plasma membrane.
	The simulated lifetimes for these processes was recorded based upon the co-diffusion category. For co-diffusion category 1, the time spent in states 1 and 2 before state 3 was measured as process 1, the time spent in state 3 was measured as process 2, and the time spent in states 1 and 2 after state 3 was measured as process 3.
	For category 2, the time spent in states 1 and 2 was measured as process 1 and the time spent in state 3 was measured as process 4.
	For category 3, the time spent in state 3 was measured as process 2 and the time spent in states 1 and 2 was measured as process 3.
	For category 4, the time spent in state 3 was measured as process 4.
	Since the activation rate (process 1) depends on the proportion of simulations beginning in state 1 (f1) versus state 2 (f2), the distributions for process 1 beginning in state 1 and process 1 beginning in state 2 were weighted (f1/(f1 + f2) and f2/(f1+f2) respectively) to create the final distribution for process 1 times.

##Score function:
	The experimental and simulated lifetime distributions were binned between 230 and 1430 ms using 60 ms bin widths. The experimental and simulated process time distributions were binned between 0 ms and 1500 ms with 20 ms bin widths.
	The score function consisted of the mean-squared error (MSE) between the experimental and simulated lifetime distributions (3), process time distributions (4), and the  co-diffusion category breakdowns. These eight values were summed for a final error score.

##Evolutionary algorithm:
	The model was parameterized using an evolutionary algorithm (EA) with 29 independent runs of 100 individuals over 40 generations. The mutation rate was set to 0.1 and the crossover rate was set to 0.5. The rate constants k1, k2, … k8 were given an allowable range between 10-4 and 20 s-1. The proportion beginning in the open tagSrc configuration with Binder was experimentally measured to be 24% (f3 = 0.24), thus the proportion of simulations beginning in the closed tagSrc configuration (f1) was given an allowable range between 10-4 and 0.76, and f2 (proportion beginning in the open tagSrc configuration) was set to 1 – f1 – f3. The EA aimed to minimize the total error as described above.

## Independence of results from tagSrc or Binder expression (`show_dox_indept`)
Loads in the process time data for the 0pg, 10pg, 100pg, and 250pg doxycycline-treated co-diffusion events.
