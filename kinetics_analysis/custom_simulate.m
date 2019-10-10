function [i,hadBinderFor12Frames,category,X] = custom_simulate(mc,numSteps,varargin)
% 20190606, Mike Pablo:
% custom_simulate is a modified version of simulate.
% The original documentation for simulate is below.
% ---------------------------------------------------------
% SIMULATE Simulate Markov chain state walks
%
% Syntax:
%
%   X = simulate(mc,numSteps)
%   X = simulate(...,param,val,...)
%
% Description:
%
%   SIMULATE computes data X on random walks through sequences of states in
%   discrete-time Markov chain mc.
%
% Input Arguments:
%
%   mc - Discrete-time Markov chain (@dtmc) object, with NumStates states.
%
%   numSteps - Positive integer specifying the number of discrete time
%        steps in each simulation.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME	VALUE
%
%   'X0'	Initial states of simulations. X0 is a nonnegative vector of
%           length NumStates, giving counts for the number of simulations
%           to begin in each state. The default is a single simulation
%           beginning from a random initial state.
%
% Output Argument:
%
%   X - X is an array of state numbers of size (1+numSteps)-by-numSims,
%       where numSims = sum(X0). The first row contains initial states.
%       Columns, in order, are all simulations beginning in the first
%       state, then all simulations beginning in the second state, etc.
%
% Notes:
%
%   o To start N simulations from state k, use:
%
%     X0 = zeros(1,NumStates);
%     X0(k) = N;
%
%   o Use DTMC/SIMPLOT to visualize the data created by SIMULATE.
%
% Example:
%
%   mc = mcmix(4,'Zeros',8);
% 
%   % Run 3 simulations of length 10 from each state:
% 
%   X = simulate(mc,10,'X0',3*ones(1,4));
%   simplot(mc,X)
%
% See also DTMC/SIMPLOT, DTMC/REDISTRIBUTE

% Copyright 2017 The MathWorks, Inc.

P = mc.P;
numStates = mc.NumStates;

% Parse inputs and set defaults:

parseObj = inputParser;

addRequired(parseObj,'numSteps',...
        @(x)validateattributes(x,{'numeric'},...
        {'scalar','positive','integer'}))
   
addParameter(parseObj,'X0',[],...
        @(x)validateattributes(x,{'numeric'},...
        {'vector','nonnegative','integer','numel',numStates}))
   
addParameter(parseObj,'dt',...
        @(x)validateattributes(x,{'numeric'},...
        {'scalar','positive'}))
    
parse(parseObj,numSteps,varargin{:});

numSteps = parseObj.Results.numSteps;
X0 = parseObj.Results.X0;
dt = parseObj.Results.dt;

if isempty(X0) % Pick random initial state
    
    X0 = zeros(1,numStates);
    p = randperm(numStates,1);
    X0(p) = 1;
    
end

numSims = sum(X0);
    
X = zeros(1+numSteps,numSims);
hadBinderFor12Frames = false;
consec_binder_time = 0;
max_consec_binder_time = 0;

startTogether = false; % if not starting tegoether, tagSrc first
endTogether = false; % if not ending together, Binder first
category=nan;

mcdt_per_frame = round(0.02/dt);

assert(0.02/dt == round(0.02/dt),'markov chain dt is non-divisible by 0.02 s');

if X0(3)
    consec_binder_time = consec_binder_time + 1;
    max_consec_binder_time = max(max_consec_binder_time,consec_binder_time);
    
    startTogether = true;
end

for j = 1:numSims
    
    % Initalize 
    simState = find(X0~=0,1);
    X0(simState) = X0(simState)-1;
    X(1,j) = simState;
    
    for i = 2:(1+numSteps)
        % Determine next state 
        u = rand;
        last_state = simState;
        simState = find(u < cumsum(P(simState,:)),1);
        X(i,j) = simState;
        
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
            else
                error('??')
            end
            break;
        end
    end
end
if i==(1+numSteps)
    i=Inf;
end
    