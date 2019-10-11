function [i,hadBinderFor12Frames,category,X] = markov_sim(P,numSteps,X0,dt)
% 20191010, Mike Pablo:
% Inputs:
%   P : Transition probability matrix (nstate x nstate)
%   numsteps: # simulation steps to perform
%   X0: one-hot vector (nstatex1) denoting what state to start from
%   dt: simulation dt.
% Outputs:
%   i: the simulation step index reached upon quitting out.
%   hadBinderFor12Frames: (true/false), was the Binder present for >12 experimental frames (240 ms)
%   category: {1,2,3,4}, co-diffusion sequence category.
%   X: {s1,s2,...,sj}, temporal sequence of Markov states until absorption.

numStates = size(P,1);

X = zeros(1+numSteps,1);
hadBinderFor12Frames = false;
consec_binder_time = 0;
max_consec_binder_time = 0;

startTogether = false; % if not starting together, tagSrc first
endTogether = false; % if not ending together, Binder first
category=nan;

mcdt_per_frame = round(0.02/dt);

assert(0.02/dt == round(0.02/dt),'markov chain dt is non-divisible by 0.02 s');

if X0(3)
    consec_binder_time = consec_binder_time + 1;
    max_consec_binder_time = max(max_consec_binder_time,consec_binder_time);

    startTogether = true;
end


% Initalize
simState = find(X0,1);
X(1) = simState;

for i = 2:(1+numSteps)
    % Determine next state
    u = rand;
    last_state = simState;
    simState = find(u < cumsum(P(simState,:)),1);
    X(i) = simState;

    if simState == 3
        consec_binder_time = consec_binder_time + 1;
        max_consec_binder_time = max(max_consec_binder_time,consec_binder_time);
    else
        consec_binder_time = 0;
    end

    if simState == 4
        if max_consec_binder_time >= 12*mcdt_per_frame
            hadBinderFor12Frames = true;
        end
        if (last_state == 3)
            endTogether = true;
        else
            endTogether = false;
        end

        % Codiffusion categories
        if ~startTogether && ~endTogether
            category = 1;
        elseif ~startTogether && endTogether
            category = 2;
        elseif startTogether && ~endTogether
            category = 3;
        elseif startTogether && endTogether
            category = 4;
        end
        break;
    end
end
if i==(1+numSteps)
    i=Inf;
end

end
