% CLUSTER_REFINEMENT_I
% Searches through tracks not assigned to clusters in the initial Voronoi
% segmentation step, then assigns those tracks to clusters if they diffuse
% into a cluster's Voronoi polygon.
%
% If a track happens to be found within more than one Voronoi polygon,
% make a warning and but assign it to the cluster it spends more time in.
% Random assignment in case of a tie.
%
% Transforms the Voronoi polygons into track IDs identifying the cluster.
%
% Input:
%   datafolders: cell array of strings for all locations of raw data (*lvPALM.mat)
%   dateprefix and filesuffix: strings for file naming, e.g.
%   '20180508_'    'clus_refinement_I'
%
%   clussuffix: strings for file input, e.g.
%   'voronoi_clus'
%
% Saves two .mat files to the working directory:
%   [dateprefix 'tag' filesuffix '.mat'] and
%   [dateprefix 'bin' filesuffix '.mat']
%
% Part of the cluster_segmentation.m pipeline.
% 2018 May 8 / Mike Pablo
function cluster_refinement_I(datafolders,dateprefix,clussuffix,filesuffix)
% Load all data
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

% Load the voronoi-segmented clusters
clusfile.tag = [dateprefix 'tag' clussuffix];
clusfile.bin = [dateprefix 'bin' clussuffix];

disp('Beginning first cluster refinement step.')

fields = {'tag','bin'};
for field_id=1:2
	  cluster_track_IDs=cell(ncells,1);
    nMultVisTrackFromUnassigned = zeros(ncells,1);
    nUnassigned = zeros(ncells,1);
  	clusdat = load(clusfile.(fields{field_id}));
  	parfor i=1:ncells
  		curr_CDSfile = [all_path.(fields{field_id}){i} '/' all_file.(fields{field_id}){i}];
      curr_ROIfile = [all_path.(fields{field_id}){i} '/roi.roi'];

  		[cluster_track_IDs{i},nMultVisTrackFromUnassigned(i),nUnassigned(i)] = refine_cluster(clusdat.segmented_clusters{i},curr_CDSfile,curr_ROIfile);
  	end
  	save([dateprefix fields{field_id} filesuffix],'cluster_track_IDs','nMultVisTrackFromUnassigned','nUnassigned');
end
end

function [cluster_track_IDs,nMultVisTrackFromUnassigned,nUnassigned] = refine_cluster(vorSeg,expdatafile,ROIfile)
% also returns the total # of tracks [that were unassigned to clusters] and that visited multiple voronoi polygons,
% and the # of tracks that were unassigned to clusters.

%% Progress update
fprintf('Processing...\n\t%s\n\t%s\n',expdatafile,ROIfile);

%% Load data
expdata = load(expdatafile,'smLinked');
roidata = ReadImageJROI(ROIfile);
%% Calculate cell ROI polyshape for simulated point rejection
warning('off','all');
cell_ps = polyshape(roidata.mnCoordinates(:,1),roidata.mnCoordinates(:,2)); % For rejection
warning('on','all');

toReject = zeros(numel(expdata.smLinked),1);
centroids = zeros(numel(expdata.smLinked),2);

% Precalculate centroids, and in association, reject unusable tracks
for i=1:numel(expdata.smLinked)
	curr_cent = mean(expdata.smLinked{i}(:,1:2));
	centroids(i,:) = curr_cent;
	if ~isinterior(cell_ps,curr_cent)
		toReject(i) = true;
	else
		toReject(i) = false;
	end
end

n_clus = numel(vorSeg);
cluster_track_IDs = cell(n_clus,1);
for j=1:numel(expdata.smLinked)
    if toReject(j)
        continue;
    end
    ct = centroids(j,:)*0.1067;

    for k=1:n_clus
        % If a track centroid is inside the cluster, it was part of the original splicing step.
        if isinterior(vorSeg{k},ct)
            cluster_track_IDs{k} = [cluster_track_IDs{k},j];
            break;
        end
    end
end

% Tracks that aren't in clusters, and are within the cell ROI are examined further
diffuseToCluster_candidateIDs = false(numel(expdata.smLinked),1);
diffuseToCluster_candidateIDs(~toReject) = true; % want tracks inside cell ROI
diffuseToCluster_candidateIDs(cell2mat(cluster_track_IDs')) = false; % don't want tracks already in clusters

diffuseToCluster_candidateIDs = find(diffuseToCluster_candidateIDs); % want numerical indexing

nMultVisTrackFromUnassigned = 0;
nUnassigned = numel(diffuseToCluster_candidateIDs);
for j=1:nUnassigned
    coords = expdata.smLinked{diffuseToCluster_candidateIDs(j)}(:,1:2)*0.1067;

    cluster_assigned = false;

    % Each position is tested to see whether it visited a cluster
    cluster_visit_matrix = false(n_clus,size(coords,1));
    for k=1:n_clus
        cluster_visit_matrix(k,:) = isinterior(vorSeg{k},coords);
        if sum(cluster_visit_matrix(k,:)) > (size(coords,1)/2) % more than half the positions are inside a cluster
            cluster_track_IDs{k} = [cluster_track_IDs{k},diffuseToCluster_candidateIDs(j)];
            cluster_assigned = true;
        end
    end

    clusters_visited = logical(sum(cluster_visit_matrix,2));
    if sum(clusters_visited)>1
        nMultVisTrackFromUnassigned = nMultVisTrackFromUnassigned+1;
    end

    % If the track wasn't immediately assigned, need to determine the maximum entry along the visit matrix
    if ~cluster_assigned
        if sum(cluster_visit_matrix(:)) == 0
            continue; % no cluster assignment for this track
        else
            cluster_tallys = sum(cluster_visit_matrix,2);
            [~,favorite_cluster] = max(cluster_tallys);
            cluster_track_IDs{favorite_cluster} = [cluster_track_IDs{favorite_cluster},diffuseToCluster_candidateIDs(j)];
        end
    end
end

end
