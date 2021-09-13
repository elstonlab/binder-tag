function error_final = markov_score_simSF(k1,k2,k3,k4,k5,k6,k7,k8,f1,dt,simSFName)
% 20191001, J. Cody Herron & Mike Pablo:
% ---------------------------------------------------------
% Description:
%
%   markov_score_simSF is a score function determining the error between
%   simulated data with a known, ground-truth parameter set and a 
%   Markov chain model for Binder - TagSrc.
%
% Example:
%    markov_score_simSF(6.49,  6.63,  1.79,  2.39,0.0, 1.30,12.5, 0.17, ...
%    0.13, 0.001,'simSF123120')
%
% Input Arguments:
%    The following variables are used for the Markov chain Binder -TagSrc
%    model. Best fit values (determined by an evolutionary algorithm) are
%    given.
%
%    VAR.     DESCRIPTION                            VALUE
%    ----------------------------------------------------------------------
%    k1       Src opening                            6.49 s^-1
%    k2       Src closing                            6.63 s^-1
%    k3       tagSrc:Binder association              1.78 s^-1
%    k4       tagSrc:Binder dissociation             2.39 s^-1
%    k5       One-step Src act. + tS:B assoc.        0.0 s^-1
%    k6       One-step Src inact. + tS:B dissoc.     1.30 s^-1
%    k7       Inactive Src PM dissociation           12.5 s^-1
%    k8       Active Src PM dissociation             0.17 s^-1
%    f1       %recruited to PM, inactive             0.13 (17% of the pop.)
%
%    dt       Time step (recommended)                0.001 s
%
%    sfName   Name of sim data
%
% Experimental Data
%
%   NAME	                     DESCRIPTION
%
%   tagsrc_lifetimes.csv         Experimental tagSrc lifetimes for
%                                codiffusion scenarios (w/ Binder)
%
%   binder_lifetimes.csv         Experimental Binder lifetimes for
%                                codiffusion scenarios (w/ tagSrc)
%
%   All events histogram.csv     Experimental tagSrc lifetimes for all
%                                scenarios
%
%   event_times.mat              Measured processes of activation,
%                                inactivation, Src PM dissociation, and
%                                active Src PM dissociation
%
%   f3                           The experimentally measured proportion of
%                                PM recruitment for active tagSrc with
%                                Binder. Note that f1 and f2 must sum to
%                                1 - f3, but neither  of these could be
%                                experimentally measured.
%
% Additional scripts used
%
%   NAME	                 DESCRIPTION
%
%   RunLength.c/RunLength.m	 Author: Jan Simon, Heidelberg
%
%   custom_simulate.m        Custom dtmc code for this model
%
% Output Argument:
%
%   error_final - Returns a sum of the MSE between the measured and simulated
%       distributions of: codiffusion categories, all scenarios tagSrc
%       lifetimes, codiffusion scenarios binder lifetimes, codiffusion
%       scenarios tagSrc lifetimes, Src activation lifetimes, Src inactivation
%       lifetimes, inactive Src PM dissociation lifetimes, and active Src PM
%       dissociation lifetimes.

f3 = 0.24; %observed proportion active w/ Binder recruited to PM
f2 = 1-f3-f1; % proportion starting in each state must sum to 1

% Discrete-Time Markov Chain transition probabilities
pab = 1-exp(-k1*dt);
pac = 1-exp(-k5*dt);
pad = 1-exp(-k7*dt);
paa = 1-pab-pac-pad;

pba = 1-exp(-k2*dt);
pbc = 1-exp(-k3*dt);
pbd = 1-exp(-k8*dt);
pbb = 1-pba-pbc-pbd;

pca = 1-exp(-k6*dt);
pcb = 1-exp(-k4*dt);
pcd = 1-exp(-k8*dt);
pcc = 1-pca-pcb-pcd;

pda = 0;
pdb = 0;
pdc = 0;
pdd = 1;

% Probability matrix
P = [paa pab pac pad; ...
     pba pbb pbc pbd; ...
     pca pcb pcc pcd; ...
     pda pdb pdc pdd];

% Require all non-diagonal rates to follow Poisson assumption.
poisson_assertion = all(P(~eye(size(P)))<0.1);

% If Poission assumption not met, return large score
if poisson_assertion ~= 1
    error_final = 100;
    return;
end



% start nrealiz simulations, each from each state.
% rng(1) % use rng(1) for testing on a single rng seed
nrealiz = 5000;

fpt = zeros(3,nrealiz); % first passage time simulation time (from start to PM recruitment)
wasCodiffusion = zeros(3,nrealiz);
codifCategory = zeros(3,nrealiz);
stateVec = cell(3,nrealiz);
binder_LT = zeros(3,nrealiz);

%lifetimes for the four events/processes
sim_process_1 = zeros(3,nrealiz);
sim_process_2 = zeros(3,nrealiz);
sim_process_3 = zeros(3,nrealiz);
sim_process_4 = zeros(3,nrealiz);


addpath('RunLength_2017_04_08/');

for i=1:3
    for j=1:nrealiz

        %start in state i
        currx0 = [0 0 0 0];
        currx0(i) = 1;

        [curr_fpt,curr_wasCodif,curr_codifCat,currstateVec] = markov_sim(P,4.5/dt,currx0,dt);

        % Simplifies the currstate Vector into state order and lifetime
        [rl_b,rl_n]=RunLength(currstateVec);

        % If Binder is re-entrant or never entered, we need to throw the trajectory.
        if sum((rl_b==3))>=2
            fpt(i,j) = nan;
            wasCodiffusion(i,j) = false;
            codifCategory(i,j) = nan;
            stateVec{i,j} = nan;
            binder_LT(i,j) = nan;

        elseif sum((rl_b==3)) == 0
            fpt(i,j) = curr_fpt*dt*1000;
            wasCodiffusion(i,j) = false;
            codifCategory(i,j) = nan;
            stateVec{i,j} = currstateVec;
            binder_LT(i,j) = nan;

        else
            % time spent in state 3 is the binder lifetime
            binder_LT(i,j) = rl_n(rl_b==3)*dt*1000;

            % Measurable Processes:
            % Process 1: time until activation
            % Process 2: '' inactivation
            % Process 3: '' Src PM dissociation
            % Process 4: '' active Src PM dissociation

            % For codiffusion scenarios where tagSrc and binder both arrived
            % and left together, process 4 is the time spent in state 3
            % (binder lifetime)
            if curr_codifCat == 4
                sim_process_4(i,j) = rl_n(rl_b==3)*dt*1000;
            end

            % For codiffusion scenarios where tagSrc and binder both
            % arrived together, but binder left first: process 3 is the
            % summed time spent in states 1 and 2, process 2 is the time
            % spent in state 3
            if curr_codifCat == 3
                    sim_process_3(i,j) = sum([rl_n(rl_b==2)*dt*1000 ; rl_n(rl_b==1)*dt*1000]);
                    sim_process_2(i,j) = rl_n(rl_b==3)*dt*1000;
            end

            % For codiffusion scenarios where tagSrc arrived first and
            % binder arrived later, but both left together: process 1 is
            % the summed time spent in states 1 and 3, process 4 is the
            % time spent in state 3
            if curr_codifCat == 2
                sim_process_1(i,j) = sum([rl_n(rl_b==2)*dt*1000 ; rl_n(rl_b==1)*dt*1000]);
                sim_process_4(i,j) = rl_n(rl_b==3)*dt*1000;
            end

            %For codiffusion scenarios where tagSrc arrived first, binder
            %arrived and left, and then tagSrc left: process 1 is the sum
            %of time spent in states 1 and 2 before state 3, process 2 is
            %the time spent in state 3, and process 3 is the the sum of
            %time spent in states 1 and 2 after state 3
            if curr_codifCat == 1
                indices_1 = find(rl_b == 1);
                indices_2 = find(rl_b == 2);
                ind_3 = find(rl_b == 3);

                pre_indices_1 = indices_1(indices_1 < ind_3);
                pre_indices_2 = indices_2(indices_2 < ind_3);
                post_indices_1 = indices_1(indices_1 > ind_3);
                post_indices_2 = indices_2(indices_2 > ind_3);

                pre_1_lft = 0;
                pre_2_lft = 0;
                post_1_lft = 0;
                post_2_lft = 0;

                if ~(length(pre_indices_1) == 0)
                    pre_1_lft = sum(rl_n(pre_indices_1)*dt*1000);
                end
                if ~(length(pre_indices_2) == 0)
                    pre_2_lft = sum(rl_n(pre_indices_2)*dt*1000);
                end
                if ~(length(post_indices_1) == 0)
                    post_1_lft = sum(rl_n(post_indices_1)*dt*1000);
                end
                if ~(length(post_indices_1) == 0)
                    post_2_lft = sum(rl_n(post_indices_2)*dt*1000);
                end

                sim_process_1(i,j) = pre_1_lft+pre_2_lft;

                sim_process_2(i,j) = rl_n(rl_b==3)*dt*1000;

                sim_process_3(i,j) = post_1_lft+post_2_lft;


            end

            fpt(i,j) = curr_fpt*dt*1000;
            wasCodiffusion(i,j) = curr_wasCodif;
            codifCategory(i,j) = curr_codifCat;
            stateVec{i,j} = currstateVec;
        end
    end
end


% Only keep fpt where codiffusion occurs
wasCodiffusion=logical(wasCodiffusion);

codif_fpt_s1 = fpt(1,wasCodiffusion(1,:));
codif_fpt_s2 = fpt(2,wasCodiffusion(2,:));
codif_fpt_s3 = fpt(3,wasCodiffusion(3,:));


% Cleaning: drop out inf-duration fpt (never exited) and <240 ms fpt
codif_fpt_s1(isinf(codif_fpt_s1)) = nan;
codif_fpt_s2(isinf(codif_fpt_s2)) = nan;
codif_fpt_s3(isinf(codif_fpt_s3)) = nan;
codif_fpt_s1(codif_fpt_s1<240) = nan;
codif_fpt_s2(codif_fpt_s2<240) = nan;
codif_fpt_s3(codif_fpt_s3<240) = nan;

codif_fpt_s1(isnan(codif_fpt_s1)) = [];
codif_fpt_s2(isnan(codif_fpt_s2)) = [];
codif_fpt_s3(isnan(codif_fpt_s3)) = [];

if isempty(codif_fpt_s1) || isempty(codif_fpt_s2)  || isempty(codif_fpt_s3)
    error_final = 100;
    return;
end

codif_binder = cat(2, binder_LT(:));
codif_binder(isinf(codif_binder)) = nan;
codif_binder(codif_binder<240) = nan;
codif_binder(isnan(codif_binder)) = [];

% Process 1 (activation rate) will depend on the proportion of the
% distribution starting in state 1 versus state 2
sim_process_1_state_1 = sim_process_1(1,:);
sim_process_1_state_2 = sim_process_1(2,:);


sim_process_1_state_1(isinf(sim_process_1_state_1)) = nan;
sim_process_1_state_2(isinf(sim_process_1_state_2)) = nan;

sim_process_2(isinf(sim_process_2)) = nan;
sim_process_3(isinf(sim_process_3)) = nan;
sim_process_4(isinf(sim_process_4)) = nan;

sim_process_1_state_1(sim_process_1_state_1 == 0) = nan;
sim_process_1_state_2(sim_process_1_state_2 == 0) = nan;

sim_process_3(sim_process_3==0) = nan;

sim_process_2(sim_process_2<240) = nan;
sim_process_4(sim_process_4<240) = nan;

% Also drop out processes longer than 1500s
sim_process_1_state_1(sim_process_1_state_1 > 1500) = nan;
sim_process_1_state_2(sim_process_1_state_2 > 1500) = nan;

sim_process_2(sim_process_2>1500) = nan;
sim_process_3(sim_process_3>1500) = nan;
sim_process_4(sim_process_4>1500) = nan;

sim_process_1_state_1(isnan(sim_process_1_state_1)) = [];
sim_process_1_state_2(isnan(sim_process_1_state_2)) = [];

sim_process_2(isnan(sim_process_2)) = [];
sim_process_3(isnan(sim_process_3)) = [];
sim_process_4(isnan(sim_process_4)) = [];


[dist1,~] = histcounts(codif_fpt_s1,'normalization','probability',...
                       'binedges',230:60:1430);
[dist2,~] = histcounts(codif_fpt_s2,'normalization','probability',...
                       'binedges',230:60:1430);
[dist3,~] = histcounts(codif_fpt_s3,'normalization','probability',...
                       'binedges',230:60:1430);

[binder_codiff_dist, ~] = histcounts(codif_binder,'normalization','probability',...
                       'binedges',230:60:1430);

% True distribution adjusted for proportion starting in each state
tag_codiff_dist = dist1*f1 + dist2*f2 + dist3*f3;

% Get individual tag and binder lifetimes

% Get and clean the data for all tagSrc lifetimes
tag_lft1 = fpt(1,:);
tag_lft2 = fpt(2,:);
tag_lft3 = fpt(3,:);


tag_lft1(isinf(tag_lft1)) = nan;
tag_lft1(tag_lft1<240) = nan;
tag_lft1(isnan(tag_lft1)) = [];
tag_lft2(isinf(tag_lft2)) = nan;
tag_lft2(tag_lft2<240) = nan;
tag_lft2(isnan(tag_lft2)) = [];
tag_lft3(isinf(tag_lft3)) = nan;
tag_lft3(tag_lft3<240) = nan;
tag_lft3(isnan(tag_lft3)) = [];


[dist_tag1,~] = histcounts(tag_lft1,'normalization','probability',...
                       'binedges',230:60:1430);
[dist_tag2,~] = histcounts(tag_lft2,'normalization','probability',...
                       'binedges',230:60:1430);
[dist_tag3,~] = histcounts(tag_lft3,'normalization','probability',...
                       'binedges',230:60:1430);

dist_all_tag = dist_tag1*f1 + dist_tag2*f2 + dist_tag3*f3;


% Get the codiffusion category percentages
codiff_cat = codifCategory(:);

% Clean and exclude observations with binder lifetime < 240 ms
codiff_cat = codiff_cat(binder_LT(:) >=240);
codiff_cat = codiff_cat(~isnan(codiff_cat));


cat1 = sum(codiff_cat == 1)/length(codiff_cat);
cat2 = sum(codiff_cat == 2)/length(codiff_cat);
cat3 = sum(codiff_cat == 3)/length(codiff_cat);
cat4 = sum(codiff_cat == 4)/length(codiff_cat);

%need to adjust so that sum (1, 2) =  f1+f2 and sum(3,4) = f3
cat1_adj = cat1*(f1+f2)/(cat1+cat2);
cat2_adj = cat2*(f1+f2)/(cat1+cat2);
cat3_adj = cat3*f3/(cat3+cat4);
cat4_adj = cat4*f3/(cat3+cat4);

sim_cat = [cat1_adj,cat2_adj,cat3_adj,cat4_adj];


%Simulation Process Distributions
[sim_process_1_state_1_dist,~]=histcounts(sim_process_1_state_1,'normalization','probability','binedges',0:20:1500);
[sim_process_1_state_2_dist,~]=histcounts(sim_process_1_state_2,'normalization','probability','binedges',0:20:1500);

% Adjust activation time distribution by appropriate proportions
sim_process_1_dist = f1/(f1+f2)*sim_process_1_state_1_dist +f2/(f1+f2)*sim_process_1_state_2_dist;

[sim_process_2_dist,~]=histcounts(sim_process_2,'normalization','probability','binedges',0:20:1500);
[sim_process_3_dist,~]=histcounts(sim_process_3,'normalization','probability','binedges',0:20:1500);
[sim_process_4_dist,~]=histcounts(sim_process_4,'normalization','probability','binedges',0:20:1500);


% Get MSEs

load(simSFName); % struct var name should be simSF

mse_codiff_cat = immse(simSF.simCodiffCat, sim_cat);

mse_codiff_tag=immse(simSF.tagCodiffDist,transpose(tag_codiff_dist));
mse_codiff_binder=immse(simSF.binderCodiffDist,transpose(binder_codiff_dist));
mse_all_tag = immse(simSF.tagAllCodiffDist,transpose(dist_all_tag));

mse_process_1 = immse(simSF.simProcess1Dist, sim_process_1_dist);
mse_process_2 = immse(simSF.simProcess2Dist, sim_process_2_dist);
mse_process_3 = immse(simSF.simProcess3Dist, sim_process_3_dist);
mse_process_4 = immse(simSF.simProcess4Dist, sim_process_4_dist);

weight = 1;
error_final = double(mse_codiff_cat + mse_all_tag + mse_codiff_binder + ...
    mse_codiff_tag + mse_process_1*weight + mse_process_2*weight + ...
    mse_process_3*weight + mse_process_4*weight);

if isnan(error_final)
    error_final = 100.;
end

end
