% SIMULATE_UNIFORM_DENSITY_DISTRIBUTIONS_50X
% Simulates negative controls for the Voronoi-based cluster segmentation.
% Randomly places points within the cell ROI, corresponding to actual observed
% (and usable) trajectories.
%
% Input:
%   datafolders: cell array of strings for all locations of raw data (*lvPALM.mat)
%   dateprefix and filesuffix: strings for file naming, e.g.
%   '20180508_'    'sim_den_distr_50x'
%
% Saves two .mat files to the working directory:
%   [dateprefix 'tag' filesuffix '.mat'] and
%   [dateprefix 'bin' filesuffix '.mat']
%
% This is for thresholding purposes, and must be compared vs. the
% experimental density distributions.
%
% Part of the cluster_segmentation.m pipeline.
function simulate_uniform_density_distributions_50x(datafolders,dateprefix,filesuffix)
rng('default'); % for reproducibility

%% Get file structures (lvPALM files)
all_path.tag = {};
all_path.bin = {};
all_file.tag = {};
all_file.bin = {};
for i=1:numel(datafolders)
    [curr_path,curr_file] = batchGetPath(datafolders{i},'RightCali_lvPALM','mat');
    all_path.tag = [all_path.tag;curr_path];
    all_file.tag = [all_file.tag;curr_file];
    [curr_path,curr_file] = batchGetPath(datafolders{i},'Left_lvPALM','mat');
    all_path.bin = [all_path.bin;curr_path];
    all_file.bin = [all_file.bin;curr_file];
end
assert(numel(all_path.tag)==numel(all_path.bin));
ncells = numel(all_path.tag);

nsimreps=50;
disp('Simulating uniform density distributions')
fields = {'tag','bin'};
for field_id=1:2
    density_distrs = cell(ncells,1);
    pad_size = 0.1; % um; for ignoring voronoi polygons near the edge.
                    % edge polygon sizes can be artificially small.
    parfor i=1:ncells
       curr_expdatafile = [all_path.(fields{field_id}){i} '/' all_file.(fields{field_id}){i}];
       curr_ROIfile = [all_path.(fields{field_id}){i} '/roi.roi'];

       density_distrs{i} = calculate_density_distributions(curr_expdatafile,curr_ROIfile,pad_size,nsimreps);
    end
  save([dateprefix fields{field_id} filesuffix],'density_distrs','pad_size')
end


end

function density_distr = calculate_density_distributions(expdatafile,roidatafile,pad_size,nsims)
% Uses Voronoi-based approach to estimate local track density.
%	density_distr is a vector is the density for each Voronoi polygon,
%       in localizations/nm2
% We reject polygons from analysis if their centroids are within 2*pad_size
%   of the cell edge, since they might have artificially high density.
% Returns 50x simulations of the density distributions.
disp(['Calculating simulated density distributions for ',expdatafile]);
expdata = load(expdatafile);
roidata = ReadImageJROI(roidatafile);
warning('off','all');
cell_ps = polyshape(roidata.mnCoordinates(:,1)*0.1067,roidata.mnCoordinates(:,2)*0.1067); % For rejection
warning('on','all');
centroids = nan(numel(expdata.smLinked),2); % x,y
window_min = min(cell_ps.Vertices);
window_max = max(cell_ps.Vertices);

density_distr = cell(nsims,1);

for sim_rep_idx = 1:nsims
    for i=1:numel(expdata.smLinked)
        % Rejection sampling: place a point randomly for every track.
        % If outside of the cell boundaries, try again.
        currcentroid = [unifrnd(window_min(1),window_max(1)),unifrnd(window_min(2),window_max(2))];
        while ~isinterior(cell_ps,currcentroid(1),currcentroid(2))
            currcentroid = [unifrnd(window_min(1),window_max(1)),unifrnd(window_min(2),window_max(2))];
        end
        centroids(i,:) = currcentroid;
    end

    % Calculate Voronoi cells and local densities, constrained by cell ROI; um units
    [v,c,~]=VoronoiLimit(centroids(:,1),centroids(:,2),...
                         'bs_ext',roidata.mnCoordinates*0.1067,...
                         'figure','off');
    npolys = size(c,1);
    xyden = nan(npolys,3);
    edgetol = 2*pad_size; % um; Exclude voronoi cell centroids near edge
                          %     they can cause artifacts in the analysis
    for j=1:npolys
        warning('off','all');
        ps = polyshape(v(c{j},:));
        warning('on','all');
       [centx,centy] = centroid(ps);
       cent_to_edge = abs(p_poly_dist(centx,centy,...
                          roidata.mnCoordinates(:,1)*0.1067,...
                          roidata.mnCoordinates(:,2)*0.1067));
       if cent_to_edge < edgetol
           continue;
       else
           den = 1/area(ps);
           % Density is now in localizations/um2
           xyden(j,:) = [centx,centy,den];
       end
    end
    to_keep = ~isnan(xyden(:,1));
    xyden = xyden(to_keep,:);

    xyden(:,3) = xyden(:,3)/(1000^2);
    density_distr{sim_rep_idx} = xyden(:,3);
end % sim reps
end
