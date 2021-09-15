% VB-HMM analysis parameter file generated by generate_runinput_files()
% Mike Pablo 2018-05-02

% Inputs
inputfile = '../vbspt_source/tag_14.mat';
trajectoryfield = 'vbspt_tracks';
% Computing strategy
parallelize_config = 0;
% Saving options
outputfile = './tag_14_HMManalysis_hidden2.mat';
jobID = 'Data from tag_14.mat :: vbspt_tracks :: 20201215';
% Data properties
timestep = 0.02;
dim = 2;
trjLmin = 2;
% Convergence and computation alternatives
runs = 25;
maxHidden = 2;
% Evaluate extra estimates including Viterbi paths
stateEstimate = 1;
maxIter = [];
relTolF = 1e-8;
tolPar = [];
% Bootstrapping
bootstrapNum = 100;
fullBootstrap = 1;
% Limits for initial conditions
init_D = [0.0001,5];
init_tD = [0.04,0.4];
% Prior distributions
% Diffusion constants
prior_type_D = 'mean_strength';
prior_D = 1;
prior_Dstrength = 5;
% Default prior choices (according to nat. meth. 2013 paper)
prior_type_Pi = 'natmet13';
prior_piStrength = 5;
prior_type_A = 'natmet13';
prior_tD = 10*timestep;
prior_tDstrength = 2*prior_tD/timestep;
