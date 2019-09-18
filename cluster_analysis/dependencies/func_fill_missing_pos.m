%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in missing coordinates of single particle tracking data caused by
% blinking. Simply insert NaNs or do interpolations
% Input:  xyt: [x1 y1 t1; x2 y2 t2: ...]
%         methods: 'nan', 'interpolation'
% Note: 
%      1. each row of xyt is not limited by three elements. [X, Y, ...., T]
%      2. Interpolation method: https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans

% Note from Bei Liu (11/16/2017); 
% A typical usage of this function for the .allfeature field in the SPT
% data would be:
%
% allfeature = trackStruct{1}.allfeature;
% allfeatureNew = func_fill_missing_pos(allfeature, 1); % interpolate missing positions
% filled_positions = allfeatureNew(:, 1:2);
%

function xyt_new = func_fill_missing_pos(xyt, varargin)

if nargin == 1
    doInterpolation = 0;
elseif nargin  == 2
    doInterpolation  = varargin{1};
else
    error('Maximum of two variable are accepted!!');
end
if (xyt(end, end) - xyt(1, end) + 1) ~=  size(xyt, 1)
    xyt_new = nan(xyt(end, end) - xyt(1, end) + 1, size(xyt, 2));
    newtt = (xyt(1, end):xyt(end, end))';
    [~, Locb] = ismember(xyt(:, end), newtt);
    xyt_new(Locb, :) = xyt;
    if doInterpolation
        xyt_new(:, 1:2) = inpaint_nans(xyt_new(:, 1:2), 2);
    end
else
    xyt_new = xyt;
end

