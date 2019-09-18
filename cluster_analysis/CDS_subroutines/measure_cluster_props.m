% MEASURE_CLUSTER_PROPS
% Measures cluster properties, including various estimates of cluster size,
% cluster region, # tracks per cluster, slow and fast diffusive zones...
%
% Input:
%   datafolders: cell array of strings for all locations of raw data (*lvPALM.mat)
%   dateprefix and filesuffix: strings for file naming, e.g.
%   '20180508_'    'clus_refinement_I'
%
%   refinementIIsuffix: strings for file input, e.g.
%   'clus_refinement_II'
%
%   vbSPTloc, vbSPTmetaloc: folders for vbSPT result and metadata, e.g.
%
%   SPTacqTime: time in seconds of SPT acquisition. Used for setting cluster lifetime to Inf.
%
% Saves a .mat files to the working directory:
%   [dateprefix filesuffix '.mat']
%
% Part of the cluster_segmentation.m pipeline.
% 2018 May 9 / Mike Pablo


function measure_cluster_props(dataloc,vbSPTloc,vbSPTmetaloc,dateprefix,suffix_refinementII_clusters,filesuffix)

[all_path,all_file] = get_file_structures(dataloc);
ncells = numel(all_path.tag);

clusID = load([dateprefix suffix_refinementII_clusters]);

tag_prop = cell(ncells,1);
bin_prop = cell(ncells,1);

[vbSPTfile,vbSPTmetafile] = get_vbSPT_file_structures(ncells,vbSPTloc,vbSPTmetaloc);

SPTacqTime = 20;
warning('Hardcoded SPT acquisition time as 20 seconds for adhesion dataset.')

parfor i=1:ncells
    curr_expdatafile_tag = [all_path.tag{i} '/' all_file.tag{i}];
    curr_expdatafile_bin = [all_path.bin{i} '/' all_file.bin{i}];

    curr_vbSPTfile_tag = vbSPTfile.tag{i};
    curr_vbSPTfile_bin = vbSPTfile.bin{i};
    curr_vbSPTmetafile_tag = vbSPTmetafile.tag{i};
    curr_vbSPTmetafile_bin = vbSPTmetafile.bin{i};

    curr_ROIfile = [all_path.tag{i} '/roi.roi'];

    tag_prop{i} = calculate_cluster_properties(clusID.tag_refined_cluster_track_IDs{i},...
                                               curr_expdatafile_tag,curr_ROIfile,...
                                               curr_vbSPTfile_tag,curr_vbSPTmetafile_tag,SPTacqTime);

    bin_prop{i} = calculate_cluster_properties(clusID.bin_refined_cluster_track_IDs{i},...
                                               curr_expdatafile_bin,curr_ROIfile,...
                                               curr_vbSPTfile_bin,curr_vbSPTmetafile_bin,SPTacqTime);
end
save([dateprefix filesuffix],'tag_prop','bin_prop');


end

function [vbSPTfile,vbSPTmetafile] = get_vbSPT_file_structures(ncells,vbSPTdata_dir,vbSPTmeta_dir)
vbSPTfile.tag = cell(ncells,1);
vbSPTfile.bin = cell(ncells,1);
vbSPTmetafile.tag = cell(ncells,1);
vbSPTmetafile.bin = cell(ncells,1);

for i=1:ncells % [assumes the vbSPT data and metadata are structured in a standarized format...]
    vbSPTfile.tag{i} = sprintf('%stag_%02d_HMManalysis_hidden2.mat',vbSPTdata_dir,i);
    vbSPTfile.bin{i} = sprintf('%sbin_%02d_HMManalysis_hidden2.mat',vbSPTdata_dir,i);
    vbSPTmetafile.tag{i} = sprintf('%stag_%02d.mat',vbSPTmeta_dir,i);
    vbSPTmetafile.bin{i} = sprintf('%sbin_%02d.mat',vbSPTmeta_dir,i);
end


end

function properties = calculate_cluster_properties(clusterTrackIDs,expdatafile,ROIfile,vbSPTfile,vbSPTlinkagefile,SPTacqTime)
%% Progress update
fprintf('Calculating cluster properties for...\n\t%s\n\t%s\n',expdatafile,ROIfile);

%% Load data
expdata = load(expdatafile);
roidata = ReadImageJROI(ROIfile);
%% Calculate cell ROI polyshape for simulated point rejection
warning('off','all');
cell_ps = polyshape(roidata.mnCoordinates(:,1),roidata.mnCoordinates(:,2)); % For rejection
warning('on','all');
numads = numel(clusterTrackIDs);

n_all_tracks = numel(expdata.corrDataStruct);
expdata.vbspt = cell(n_all_tracks,1);
vbsptdata = load(vbSPTfile,'Wbest');
vbsptlinkagedata = load(vbSPTlinkagefile,'vbspt_trackIDs');

SPTacqFrame = round(SPTacqTime/0.02); % (20 ms per frame)

% Annotate the expdata structure with its corresponding vbSPT entries
n_vbsptTracks = numel(vbsptlinkagedata.vbspt_trackIDs);
for i=1:n_vbsptTracks
    trackID = vbsptlinkagedata.vbspt_trackIDs(i);
    non_interp_vbspt = vbsptdata.Wbest.est2.sMaxP{i};
    observed_frames = expdata.corrDataStruct{trackID}.allfeature(:,end) - ...
                      expdata.corrDataStruct{trackID}.startFrame + 1;
    % Final frame does not have a diffusive state.
    expdata.vbspt{trackID} = non_interp_vbspt(observed_frames(1:end-1));
end

toReject = zeros(n_all_tracks,1);
centroids = zeros(n_all_tracks,2);

% Precalculate centroids, and in association, reject unusable tracks
for i=1:n_all_tracks
	curr_cent = mean(expdata.corrDataStruct{i}.position);
	centroids(i,:) = curr_cent;
	if ~isinterior(cell_ps,curr_cent)
		toReject(i) = true;
	else
		toReject(i) = false;
	end
end

linearized_clusterTrackIDs =[clusterTrackIDs{:}];
n_clus = numel(linearized_clusterTrackIDs);
n_clus_per_adFrame = cellfun(@numel,clusterTrackIDs);

start_clus_ID_per_adFrame = [1;cumsum(n_clus_per_adFrame(1:end-1))+1];
end_clus_ID_per_adFrame = [cumsum(n_clus_per_adFrame(1:end))];

clus_indexing_per_adFrame = [start_clus_ID_per_adFrame,end_clus_ID_per_adFrame];

properties.trackIDs = cell(n_clus,1);
properties.N = zeros(n_clus,1);

properties.cent.fitR = zeros(n_clus,1); % fit enclosing circle to cluster centroid
properties.cent.COM = zeros(n_clus,2);   % track centroid-based center of mass

properties.allLoc.KDEeffR50 = zeros(n_clus,1); % effective radius from all regions at 50% max via kernel density estimation
properties.allLoc.KDEeffR95 = zeros(n_clus,1); % effective radius from all regions at 5% max via kernel density estimation
properties.allLoc.KDEcirc50 = zeros(n_clus,1); % circularity from all regions at 50% max via kernel density estimation
properties.allLoc.KDEcirc95 = zeros(n_clus,1); % circularity from all regions at 5% max via kernel density estimation
properties.allLoc.KDEeccen50 = zeros(n_clus,1); % eccentricity from all regions at 50% max via kernel density estimation
properties.allLoc.KDEeccen95 = zeros(n_clus,1); % eccentricity from all regions at 5% max via kernel density estimation
properties.allLoc.KDEarea50 = zeros(n_clus,1); % area from all regions at 50% max via kernel density estimation
properties.allLoc.KDEarea95 = zeros(n_clus,1); % area from all regions at 5% max via kernel density estimation
properties.allLoc.KDEmajAx50 = zeros(n_clus,1); % majoraxis from all regions at 50% max via kernel density estimation
properties.allLoc.KDEmajAx95 = zeros(n_clus,1); % majoraxis from all regions at 5% max via kernel density estimation
properties.allLoc.KDEminAx50 = zeros(n_clus,1); % minoraxis from all regions at 50% max via kernel density estimation
properties.allLoc.KDEminAx95 = zeros(n_clus,1); % minoraxis from all regions at 5% max via kernel density estimation
properties.allLoc.KDEbound50 = cell(n_clus,1); % boundaries from 50% max via kernel density estimation
properties.allLoc.KDEbound95 = cell(n_clus,1); % boundaries from 5% max via kernel density estimation

properties.allLoc.COM = zeros(n_clus,2); % all localization-based center of mass
properties.allLoc.KDEmaskCOM50 = zeros(n_clus,2); % center of mass of the 50%KDE threshold mask
properties.allLoc.KDEmaskCOM95 = zeros(n_clus,2); % center of mass of the 95%KDE threshold mask

properties.vbSPT.slowFrac = zeros(n_clus,1);
properties.vbSPT.slowd2centCOM = cell(n_clus,1); % all distances to the cluster centroid COM.
properties.vbSPT.fastd2centCOM = cell(n_clus,1); % all distances to the cluster center.
properties.vbSPT.slowd2allLocCOM = cell(n_clus,1); % all distances to the all localization COM.
properties.vbSPT.fastd2allLocCOM = cell(n_clus,1); % all distances to the all localization center.

properties.vbSPT.sloweffR50 = zeros(n_clus,1); % effective radius at 50% max via kernel density estimation
properties.vbSPT.sloweffR95 = zeros(n_clus,1); % effective radius at 5% max via kernel density estimation
properties.vbSPT.fasteffR50 = zeros(n_clus,1); % effective radius at 50% max via kernel density estimation
properties.vbSPT.fasteffR95 = zeros(n_clus,1); % effective radius at 5% max via kernel density estimation
properties.vbSPT.slowcirc50 = zeros(n_clus,1); % circularity at 50% max via kernel density estimation
properties.vbSPT.slowcirc95 = zeros(n_clus,1); % circularity at 5% max via kernel density estimation
properties.vbSPT.fastcirc50 = zeros(n_clus,1); % circularity at 50% max via kernel density estimation
properties.vbSPT.fastcirc95 = zeros(n_clus,1); % circularity at 5% max via kernel density estimation
properties.vbSPT.sloweccen50 = zeros(n_clus,1); % eccentricity at 50% max via kernel density estimation
properties.vbSPT.sloweccen95 = zeros(n_clus,1); % eccentricity at 5% max via kernel density estimation
properties.vbSPT.fasteccen50 = zeros(n_clus,1); % eccentricity at 50% max via kernel density estimation
properties.vbSPT.fasteccen95 = zeros(n_clus,1); % eccentricity at 5% max via kernel density estimation

properties.vbSPT.slowarea50 = zeros(n_clus,1); % area at 50% max via kernel density estimation
properties.vbSPT.slowarea95 = zeros(n_clus,1); % area at 5% max via kernel density estimation
properties.vbSPT.fastarea50 = zeros(n_clus,1); % area at 50% max via kernel density estimation
properties.vbSPT.fastarea95 = zeros(n_clus,1); % area at 5% max via kernel density estimation
properties.vbSPT.slowmajAx50 = zeros(n_clus,1); % majoraxis at 50% max via kernel density estimation
properties.vbSPT.slowmajAx95 = zeros(n_clus,1); % majoraxis at 5% max via kernel density estimation
properties.vbSPT.fastmajAx50 = zeros(n_clus,1); % majoraxis at 50% max via kernel density estimation
properties.vbSPT.fastmajAx95 = zeros(n_clus,1); % majoraxis at 5% max via kernel density estimation
properties.vbSPT.slowminAx50 = zeros(n_clus,1); % minoraxis at 50% max via kernel density estimation
properties.vbSPT.slowminAx95 = zeros(n_clus,1); % minoraxis at 5% max via kernel density estimation
properties.vbSPT.fastminAx50 = zeros(n_clus,1); % minoraxis at 50% max via kernel density estimation
properties.vbSPT.fastminAx95 = zeros(n_clus,1); % minoraxis at 5% max via kernel density estimation
properties.vbSPT.slowbound50 = cell(n_clus,1); % boundary at 50% max via kernel density estimation
properties.vbSPT.slowbound95 = cell(n_clus,1); % boundary at 5% max via kernel density estimation
properties.vbSPT.fastbound50 = cell(n_clus,1); % boundary at 50% max via kernel density estimation
properties.vbSPT.fastbound95 = cell(n_clus,1); % boundary at 5% max via kernel density estimation
properties.vbSPT.slowCOM50 = zeros(n_clus,2); % mask center of mass at 50% max via kernel density estimation
properties.vbSPT.slowCOM95 = zeros(n_clus,2); % mask center of mass at 5% max via kernel density estimation
properties.vbSPT.fastCOM50 = zeros(n_clus,2); % mask center of mass at 50% max via kernel density estimation
properties.vbSPT.fastCOM95 = zeros(n_clus,2); % mask center of mass at 5% max via kernel density estimation

properties.adhesion.d2ad = zeros(n_clus,1); % distance to adhesion from cluster KDE-based measure in um
properties.adhesion.adsize = zeros(n_clus,1); % size of adhesion in px
properties.adhesion.adframe = zeros(n_clus,1); % adhesion frame ID for the cluster
properties.adhesion.adID = zeros(n_clus,1); % ID of the adhesion

properties.lifetime = zeros(n_clus,1); % observed lifetime of cluster (s)
properties.rec_int = cell(n_clus,1); % recruitment intervals (s)

for i=1:numads
    curr_clus_start = clus_indexing_per_adFrame(i,1);
    curr_clus_end = clus_indexing_per_adFrame(i,2);
    properties.adhesion.adframe(curr_clus_start:curr_clus_end) = i;
end
for i=1:n_clus
   properties.trackIDs{i} = linearized_clusterTrackIDs{i};
   properties.N(i) = numel(linearized_clusterTrackIDs{i});
   clus_track_centroids = centroids(linearized_clusterTrackIDs{i},:)*0.1067;
   allLocalizations = cell2mat(cellfun(@(x) x.position,expdata.corrDataStruct(linearized_clusterTrackIDs{i}),'uniformoutput',false))*0.1067;

   properties.cent.COM(i,:) = mean(clus_track_centroids);
   properties.allLoc.COM(i,:) = mean(allLocalizations);

   % Centroid-based circle fit.
   if size(clus_track_centroids,1) == 1
       % zero radius
       properties.cent.fitR(i) = 0;
   elseif size(clus_track_centroids,1) == 2
       % half distance between points
       properties.cent.fitR(i) = dist_to_point(clus_track_centroids(1,:),...
                                               clus_track_centroids(2,:))/2;
   elseif size(clus_track_centroids,1) >= 3
       % fit minimum enclosing circle
       [centFitCircR,~,~] = ExactMinBoundCircle(clus_track_centroids);
       properties.cent.fitR(i) = centFitCircR;
   end

   % All localization KDE radius fits [KDE domain is on a (0.01 um)^2 grid]
   [area50,r50,circ50,majAx50,minAx50,eccen50,bound50,com50,...
    area95,r95,circ95,majAx95,minAx95,eccen95,bound95,com95] = fit_KDE_radius(allLocalizations);
   properties.allLoc.KDEeffR50(i) = r50;
   properties.allLoc.KDEeffR95(i) = r95;
   properties.allLoc.KDEcirc50(i) = circ50;
   properties.allLoc.KDEcirc95(i) = circ95;
   properties.allLoc.KDEeccen50(i) = eccen50;
   properties.allLoc.KDEeccen95(i) = eccen95;
   properties.allLoc.KDEarea50(i) = area50;
   properties.allLoc.KDEarea95(i) = area95;
   properties.allLoc.KDEmajAx50(i) = majAx50;
   properties.allLoc.KDEmajAx95(i) = majAx95;
   properties.allLoc.KDEminAx50(i) = minAx50;
   properties.allLoc.KDEminAx95(i) = minAx95;
   properties.allLoc.KDEbound50{i} = bound50;
   properties.allLoc.KDEbound95{i} = bound95;
   properties.allLoc.KDEmaskCOM50(i,:) = com50;
   properties.allLoc.KDEmaskCOM95(i,:) = com95;

   % vbSPT analyses
   all_states = cell2mat(cellfun(@(x) x,expdata.vbspt(properties.trackIDs{i}),'uniformoutput',false));
   properties.vbSPT.slowFrac(i) = sum(all_states==1)./numel(all_states);

   all_slow_localizations = cell2mat(cellfun(@(x,y) x.position(y==1,1:2)*0.1067,expdata.corrDataStruct(properties.trackIDs{i}),expdata.vbspt(properties.trackIDs{i}),'uniformoutput',false));
   all_fast_localizations = cell2mat(cellfun(@(x,y) x.position(y==2,1:2)*0.1067,expdata.corrDataStruct(properties.trackIDs{i}),expdata.vbspt(properties.trackIDs{i}),'uniformoutput',false));

   % distance to track centroid-based cluster COM
   properties.vbSPT.slowd2centCOM{i} = dist_to_point(all_slow_localizations,...
                                                    properties.cent.COM(i,:));
   properties.vbSPT.fastd2centCOM{i} = dist_to_point(all_fast_localizations,...
                                                    properties.cent.COM(i,:));

   % distance to all track localization-based cluster COM
   properties.vbSPT.slowd2allLocCOM{i} = dist_to_point(all_slow_localizations,...
                                                    properties.allLoc.COM(i,:));
   properties.vbSPT.fastd2allLocCOM{i} = dist_to_point(all_fast_localizations,...
                                                    properties.allLoc.COM(i,:));

   % kernel density estimate estimates of slow and fast zone sizes
   [area50,r50,circ50,majAx50,minAx50,eccen50,bound50,com50,...
    area95,r95,circ95,majAx95,minAx95,eccen95,bound95,com95] = fit_KDE_radius(all_slow_localizations);
   properties.vbSPT.sloweffR50(i) = r50;
   properties.vbSPT.sloweffR95(i) = r95;
   properties.vbSPT.slowcirc50(i) = circ50;
   properties.vbSPT.slowcirc95(i) = circ95;
   properties.vbSPT.sloweccen50(i) = eccen50;
   properties.vbSPT.sloweccen95(i) = eccen95;
   properties.vbSPT.slowarea50(i) = area50;
   properties.vbSPT.slowarea95(i) = area95;
   properties.vbSPT.slowmajAx50(i) = majAx50;
   properties.vbSPT.slowmajAx95(i) = majAx95;
   properties.vbSPT.slowminAx50(i) = minAx50;
   properties.vbSPT.slowminAx95(i) = minAx95;
   properties.vbSPT.slowbound50{i} = bound50;
   properties.vbSPT.slowbound95{i} = bound95;
   properties.vbSPT.slowCOM50(i,:) = com50;
   properties.vbSPT.slowCOM95(i,:) = com95;

   [area50,r50,circ50,majAx50,minAx50,eccen50,bound50,com50,...
    area95,r95,circ95,majAx95,minAx95,eccen95,bound95,com95] = fit_KDE_radius(all_fast_localizations);
   properties.vbSPT.fasteffR50(i) = r50;
   properties.vbSPT.fasteffR95(i) = r95;
   properties.vbSPT.fastcirc50(i) = circ50;
   properties.vbSPT.fastcirc95(i) = circ95;
   properties.vbSPT.fasteccen50(i) = eccen50;
   properties.vbSPT.fasteccen95(i) = eccen95;
   properties.vbSPT.fastarea50(i) = area50;
   properties.vbSPT.fastarea95(i) = area95;
   properties.vbSPT.fastmajAx50(i) = majAx50;
   properties.vbSPT.fastmajAx95(i) = majAx95;
   properties.vbSPT.fastminAx50(i) = minAx50;
   properties.vbSPT.fastminAx95(i) = minAx95;
   properties.vbSPT.fastbound50{i} = bound50;
   properties.vbSPT.fastbound95{i} = bound95;
   properties.vbSPT.fastCOM50(i,:) = com50;
   properties.vbSPT.fastCOM95(i,:) = com95;

   % lifetimes and recruitment intervals
   track_starts = cellfun(@(x) x.startFrame,expdata.corrDataStruct(linearized_clusterTrackIDs{i}));
   track_ends   = cellfun(@(x) x.endFrame,expdata.corrDataStruct(linearized_clusterTrackIDs{i}));
   track_starts = sort(track_starts);
   properties.rec_int{i} = (track_starts(2:end)-track_starts(1:end-1))*0.02;
   if mod(max(track_ends),SPTacqFrame) == 0
       properties.lifetime(i) = Inf; % cannot measure lifetime
   else
       properties.lifetime(i) = (max(track_ends) - min(track_starts))*0.02;
   end

   % adhesion assignment
   curr_clus_COM = properties.allLoc.KDEmaskCOM50(i,:);
   curr_adFrame = properties.adhesion.adframe(i);
   curr_numads = numel(expdata.focalAdhesionAll{curr_adFrame});
   curr_dist2ads = zeros(curr_numads,1);
   curr_adhesionID = nan;
   for j=1:curr_numads
       focalAdCoords = expdata.focalAdhesionAll{curr_adFrame}{j}.sROI.mnCoordinates*0.1067;
       curr_dist2ads(j) = p_poly_dist(curr_clus_COM(1),curr_clus_COM(2),...
                                      focalAdCoords(:,1),focalAdCoords(:,2));
       if curr_dist2ads(j) < 0
           curr_adhesionID = j;
           minAdDist = curr_dist2ads(j);
           break;
       end
   end
   if isnan(curr_adhesionID)
       [minAdDist,curr_adhesionID] = min(curr_dist2ads);
   end
   properties.adhesion.adID(i) = curr_adhesionID;
   properties.adhesion.d2ad(i) = minAdDist;
   properties.adsize(i) = numel(expdata.focalAdhesionAll{curr_adFrame}{curr_adhesionID}.pixelListInd);
end
end


function dist = dist_to_point(vec,point)
if isempty(vec) || isempty(point)
    dist = nan;
else
    dist = sqrt( (vec(:,1)-point(1)).^2 + (vec(:,2)-point(2)).^2 );
end
end

function [area50,r50,circ50,majAx50,minAx50,eccen50,bound50,com50,...
          area95,r95,circ95,majAx95,minAx95,eccen95,bound95,com95] = fit_KDE_radius(points)
% performs a 2D kernel density estimate to points, a Nx2 matrix of coordinates in um..
% thresholds for 50% max and 5% max, and
% then returns areas, effective radii, circularities, major axes, minor axes, eccentricities, and
% perimeter coordinates corresponding to those regions.
%
% if more than one region is detected after thresholding, returns nan values.
%
% the KDE is estimated on a region padded by 1 um from observed extremes, on a 0.01 um spacing
if isempty(points) || (size(points,1) == 1) % no points, or only one entry
    r50=nan;
    r95=nan;
    circ50 = nan;
    circ95 = nan;
    eccen50 = nan;
    eccen95 = nan;
    area50 = nan;
    area95 = nan;
    majAx50=nan;
    majAx95 = nan;
    minAx50=nan;
    minAx95=nan;
    bound50=[nan nan];
    bound95=[nan nan];
    com50=[nan nan];
    com95=[nan nan];
else
    gridspacing = 0.01; % um/px
    domain = [min(points(:,1))-1, max(points(:,1))+1, ...
              min(points(:,2))-1, max(points(:,2))+1];

    gridx = domain(1):gridspacing:domain(2);
    gridy = domain(3):gridspacing:domain(4);
    [mx,my]=meshgrid(gridx,gridy); % mesh representation
    vx = mx(:); % vector representation
    vy = my(:);
    vxy = [vx vy];

    vks = ksdensity(points,vxy);
    maxvks = max(vks(:));
    bin50 = vks > (maxvks*0.50);
    bin95 = vks > (maxvks*0.05);

    um2_per_px = gridspacing^2;

    % calculate area first in px for the circularity measurement
    area50 = sum(bin50(:));
    area95 = sum(bin95(:));

    % calculate circularity via perimeter; set a nan value if disjointed region
    % but keep the effective radius measurement from above.
    imbin50 = reshape(bin50,length(gridy),length(gridx));
    imbin95 = reshape(bin95,length(gridy),length(gridx));
    % note --> display would be imagesc(grix,gridy,imbin50). Recall that y = rows, x = columns.

    perim50 = regionprops(imbin50,'perimeter');
    perim95 = regionprops(imbin95,'perimeter');

    if numel(perim50) > 1 || numel(perim50) == 0
        circ50 = nan;
        eccen50 = nan;
        majAx50 = nan;
        minAx50 = nan;
        com50 = [nan nan];
        bound50 = [nan nan];
        %warning('non-continuous region in 50% measure')
    else
        circ50 = 4*pi*area50/(perim50.Perimeter)^2;
        eccen50 = regionprops(imbin50,'eccentricity');
        eccen50 = eccen50.Eccentricity;

        majAx50 = regionprops(imbin50,'MajorAxisLength');
        minAx50 = regionprops(imbin50,'MinorAxisLength');
        majAx50 = majAx50.MajorAxisLength*gridspacing; % convert to microns
        minAx50 = minAx50.MinorAxisLength*gridspacing;

        com50 = regionprops(imbin50,'Centroid');

        com50 = com50.Centroid;
        % convert local px-based center of mass to global coordinate system in micron
        com50(1) = (com50(1)-1)*gridspacing+gridx(1);
        com50(2) = (com50(2)-1)*gridspacing+gridy(1);

        % get the exterior boundary, and convert these to global coordinates
        B50=bwboundaries(imbin50,'noholes');
        bound50 = [gridx(B50{1}(:,2))',gridy(B50{1}(:,1))'];
    end
    if numel(perim95) > 1 || numel(perim95) == 0
        circ95 = nan;
        eccen95 = nan;
        majAx95 = nan;
        minAx95 = nan;
        com95 = [nan nan];
        bound95 = [nan nan];
        %warning('non-continuous region in 95% measure')
    else
        circ95 = 4*pi*area95/(perim95.Perimeter)^2;
        eccen95 = regionprops(imbin95,'eccentricity');
        eccen95 = eccen95.Eccentricity;

        majAx95 = regionprops(imbin95,'MajorAxisLength');
        minAx95 = regionprops(imbin95,'MinorAxisLength');
        majAx95 = majAx95.MajorAxisLength*gridspacing;
        minAx95 = minAx95.MinorAxisLength*gridspacing;

        com95 = regionprops(imbin95,'Centroid');
        com95 = com95.Centroid;
        com95(1) = (com95(1)-1)*gridspacing+gridx(1);
        com95(2) = (com95(2)-1)*gridspacing+gridy(1);

        % get the exterior boundary, and convert these to global coordinates
        B95=bwboundaries(imbin95,'noholes');
        bound95 = [gridx(B95{1}(:,2))',gridy(B95{1}(:,1))'];
    end

    % convert area to squared-microns
    area50 = area50*um2_per_px;
    area95 = area95*um2_per_px;

    r50 = sqrt(area50/pi);
    r95 = sqrt(area95/pi);

	nm2_per_um2 = 1000^2;
    area50 = area50*nm2_per_um2;
    area95 = area95*nm2_per_um2;
end
end


function [all_path,all_file] = get_file_structures(datafolders)
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
end
