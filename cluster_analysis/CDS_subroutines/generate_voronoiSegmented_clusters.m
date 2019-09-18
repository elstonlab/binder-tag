% GENERATE_VORONOISEGMENTED_CLUSTERS
% Generates Voronoi segmentations from experimental data,
% then uses the 50x simulation-based threshold to delete irrelevant
% polygons. Merges remaining touching polygons to obtain clusters (polyshapes)
%
% I calculate the threshold separately for each adhesion frame.
%
% Input:
%   datafolders: cell array of strings for all locations of raw data (*lvPALM.mat)
%   dateprefix and filesuffix: strings for file naming, e.g.
%   '20180508_'    'sim_den_distr_50x'
%
%   threshsuffix: strings for file input, e.g.
%   'voronoi_thresh'
%
% Saves two .mat files to the working directory:
%   [dateprefix 'tag' filesuffix '.mat'] and
%   [dateprefix 'bin' filesuffix '.mat']
%
% Part of the cluster_segmentation.m pipeline.


function generate_voronoiSegmented_clusters(datafolders,dateprefix,filesuffix,threshsuffix)
%% Get file structures (CDS files)
all_path.tag = {};
all_path.bin = {};
all_file.tag = {};
all_file.bin = {};
for i=1:numel(datafolders)
    [curr_path,curr_file] = batchGetPath(datafolders{i},'RightCali_lvPALM_corrDataStruct','mat');
    all_path.tag = [all_path.tag;curr_path];
    all_file.tag = [all_file.tag;curr_file];
    [curr_path,curr_file] = batchGetPath(datafolders{i},'Left_lvPALM_corrDataStruct','mat');
    all_path.bin = [all_path.bin;curr_path];
    all_file.bin = [all_file.bin;curr_file];
end
assert(numel(all_path.tag)==numel(all_path.bin));
ncells = numel(all_path.tag);

segmented_clusters = cell(ncells,1);
disp('Beginning Voronoi segmentation of clusters')
fields = {'tag','bin'};
threshold_data = load([dateprefix threshsuffix]);
for field_id=1:2
	  pad_size = 0.1; % um; for ignoring voronoi polygons near the edge.
                    % edge polygon sizes can be artificially small.
    parfor i=1:ncells
       curr_CDSfile = [all_path.(fields{field_id}){i} '/' all_file.(fields{field_id}){i}];
       curr_ROIfile = [all_path.(fields{field_id}){i} '/roi.roi'];

       if strcmp(fields{field_id},'tag')
           curr_thresh = threshold_data.tag_cutoffs(i); % cell array for CDS format
       elseif strcmp(fields{field_id},'bin')
           curr_thresh = threshold_data.bin_cutoffs(i);  % cell array for CDS format
       end
       segmented_clusters{i} = segment_clusters(curr_CDSfile,curr_ROIfile,pad_size,curr_thresh);
    end
    save([dateprefix fields{field_id} filesuffix],'segmented_clusters','pad_size','threshold_data');
end


end

function all_segmented_clusters = segment_clusters(expdatafile,roidatafile,pad_size,thresholdvalue)
% Uses Voronoi-based approach to segment clusters in the expdatafile.
disp(['Performing thresholded Voronoi segmentation and cluster fusion for ',expdatafile]);
expdata = load(expdatafile);
roidata = ReadImageJROI(roidatafile);
warning('off','all');
cell_ps = polyshape(roidata.mnCoordinates(:,1)*0.1067,roidata.mnCoordinates(:,2)*0.1067); % For rejection
warning('on','all');
numads = numel(expdata.focalAdhesionAll);
all_segmented_clusters = cell(numads,1);
thresholdvalue=cell2mat(thresholdvalue);

for adframe_idx=1:numads
    curr_thresholdvalue = thresholdvalue(adframe_idx);

    current_tracks = cellfun(@(x) x.adhesionFrame==adframe_idx,expdata.corrDataStruct);
    current_tracks_numeric_idx = find(current_tracks);
    num_tracks = sum(current_tracks);
    centroids = nan(num_tracks,2); % x,y

    for i=1:num_tracks
        curr_track_idx = current_tracks_numeric_idx(i);
        currcentroid = mean(expdata.corrDataStruct{curr_track_idx}.position(:,1:2)*0.1067);
        % Some centroids might be outside of the cell boundary.
        % These most usually are associated with camera artifacts or spurious
        % signal (i.e. debris) that needed to be excluded.
        if ~isinterior(cell_ps,currcentroid(1),currcentroid(2))
            continue;
        end
        % Convert centroid to um
        centroids(i,:) = currcentroid;
    end
    % Drop out any tracks that got skipped
    centroids = centroids(~isnan(centroids(:,1)),:);

    % Calculate Voronoi cells and local densities, constrained by cell ROI; um units
    [v,c,~]=VoronoiLimit(centroids(:,1),centroids(:,2),...
                         'bs_ext',roidata.mnCoordinates*0.1067,...
                         'figure','off');

    npolys = size(c,1);
    edgetol = 2*pad_size; % um; Exclude voronoi cell centroids near edge
                          % they can cause artifacts in the analysis
    poly_to_drop = false(npolys,1);
    poly_dens = zeros(npolys,1);
    all_polys = cell(npolys,1);
    for j=1:npolys
        warning('off','all');
        ps = polyshape(v(c{j},:));
        all_polys{j} = ps;
        warning('on','all');
        [centx,centy] = centroid(ps);
        cent_to_edge = abs(p_poly_dist(centx,centy,...
                          roidata.mnCoordinates(:,1)*0.1067,...
                          roidata.mnCoordinates(:,2)*0.1067));
       if cent_to_edge < edgetol
           poly_to_drop(j) = true;
           continue;
       else
           den = 1/area(ps); % Density is now in localizations/um2
           den = den/(1000^2); % density is now in localiations/nm2
           poly_dens(j) = den;
       end
    end
    poly_to_drop = poly_to_drop | (poly_dens<curr_thresholdvalue);
    all_polys = all_polys(~poly_to_drop);
    npolys = numel(all_polys);
    % Use graph theoretic method to get connected polygons.
    adjacency_mat = false(npolys,npolys);
    for j=1:npolys
        for k=j:npolys
            v1 = all_polys{j}.Vertices;
            v2 = all_polys{k}.Vertices;
            if share_vertex(v1,v2)
                adjacency_mat(j,k) = true;
                adjacency_mat(k,j) = true;
            end
        end
    end
    % Get the IDs of base polygons composing unique clusters
    polygraph = graph(adjacency_mat);
    clusters=conncomp(polygraph); % Each polygon is assigned to a cluster ID

    cluster_list = unique(clusters);
    nclus = numel(cluster_list);
    segmented_clusters = cell(nclus,1); % we'll remove empty entries after
    for j=1:nclus
        curr_clus = cluster_list(j);
        curr_poly_in_clus = find(curr_clus==clusters);
        % Want to have at least three sufficiently dense polygons to start a cluster
        if numel(curr_poly_in_clus) < 3
            continue;
        end
        % Start merging together polygons that are attached
        curr_seg_clus = union(all_polys{curr_poly_in_clus(1)},...
                              all_polys{curr_poly_in_clus(2)});
        for k=3:numel(curr_poly_in_clus)
            curr_seg_clus = union(curr_seg_clus,...
                                  all_polys{curr_poly_in_clus(k)});
        end
        segmented_clusters{j} = curr_seg_clus;
    end
    segmented_clusters = segmented_clusters(cellfun(@(x) ~isempty(x),segmented_clusters));

    all_segmented_clusters{adframe_idx} = segmented_clusters;
end
end

function dotouch = share_vertex(v1,v2)
nv1 = size(v1,1);
dotouch=false;
for i=1:nv1
    if any(v1(i,:)==v2)
        dotouch = true;
        break;
    end
end
end
