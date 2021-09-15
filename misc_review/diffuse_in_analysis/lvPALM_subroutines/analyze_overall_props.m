% ANALYZE_OVERALL_PROPS
% Analyzes and visualizes overall properties, including % of tracks participating in clusters,
% % of in-cluster tracks that started in clusters, % of tracks that diffused into clusters,
% "overall" and "continuous" per-track cluster on-times and off-times.
%
% all analyses EXCEPT for the % tracks participating in clusters are repeated
% for each approach to establishing the cluster center and size.
% the % tracks participating is by definition based on the # tracks defining clusters.
%
% Input:
%   dateprefix and filesuffix: strings for file naming, e.g.
%   '20180508_'    'analyzed_clus_prop'
%
%   propsuffix: strings for file input, e.g.
%   'measured_cluster_props'
%
%   nMinTracks: integer minimum required # of tracks in a cluster to be included in analysis
%
% Saves a .mat file to the working directory:
%   [dateprefix filesuffix '.mat']
%
% Also produces figures for each type of analysis carried out.
%
% Part of the cluster_segmentation.m pipeline.
% 2018 May 9 / Mike Pablo

function analyze_overall_props(dataloc,vbSPTloc,vbSPTmetaloc,dateprefix,measpropsuffix,analpropsuffix,...
                               nMinTracks,filesuffix,parentFigDir)
[all_path,all_file] = get_file_structures(dataloc);
ncells = numel(all_path.tag);

measprop = load([dateprefix measpropsuffix]); % measured properties, needed for track IDs and cluster locs.
analprop = load([dateprefix analpropsuffix]); % analyzed properties, needed for cluster sizes

[vbSPTfile,vbSPTmetafile] = get_vbSPT_file_structures(ncells,vbSPTloc,vbSPTmetaloc);

overall_prop_tag = cell(ncells,1);
overall_prop_bin = cell(ncells,1);

figSubDir = 'overall_analysis/';
parentFigDir = [parentFigDir figSubDir];
warning('off');
mkdir(parentFigDir);
warning('on');

parfor i=1:ncells
    curr_expdatafile_tag = [all_path.tag{i} '/' all_file.tag{i}];
    curr_expdatafile_bin = [all_path.bin{i} '/' all_file.bin{i}];

    curr_vbSPTfile_tag = vbSPTfile.tag{i};
    curr_vbSPTfile_bin = vbSPTfile.bin{i};
    curr_vbSPTmetafile_tag = vbSPTmetafile.tag{i};
    curr_vbSPTmetafile_bin = vbSPTmetafile.bin{i};

    curr_ROIfile = [all_path.tag{i} '/roi.roi'];

    overall_prop_tag{i} = calculate_overall_properties(curr_expdatafile_tag,curr_ROIfile,...
                                               curr_vbSPTfile_tag,curr_vbSPTmetafile_tag,...
                                               measprop.tag_prop{i},analprop.propPointEst.tag{i});

    overall_prop_bin{i} = calculate_overall_properties(curr_expdatafile_bin,curr_ROIfile,...
                                               curr_vbSPTfile_bin,curr_vbSPTmetafile_bin,...
                                               measprop.bin_prop{i},analprop.propPointEst.bin{i});

end
% Struct-based storage handled outside of the parfor loop.
overall_prop.tag = overall_prop_tag;
overall_prop.bin = overall_prop_bin;
save([dateprefix filesuffix],'overall_prop');
end

function overall = calculate_overall_properties(expdatafile,ROIfile,vbSPTfile,vbSPTlinkagefile,...
                                                measprops,analprops)
fprintf('Calculating overall properties for...\n\t%s\n',expdatafile);

%% Load data
expdata = load(expdatafile,'smLinked');
roidata = ReadImageJROI(ROIfile);
%% Calculate cell ROI polyshape for simulated point rejection
warning('off','all');
cell_ps = polyshape(roidata.mnCoordinates(:,1),roidata.mnCoordinates(:,2)); % For rejection
warning('on','all');

toReject = zeros(numel(expdata.smLinked),1);
centroids = zeros(numel(expdata.smLinked),2);

expdata.vbspt = cell(numel(expdata.smLinked),1);
expdata.vbspt_interp = cell(numel(expdata.smLinked),1);
vbsptdata = load(vbSPTfile,'Wbest');
vbsptlinkagedata = load(vbSPTlinkagefile,'vbspt_trackIDs');

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

% Annotate expdata vbSPT entries
n_vbsptTracks = numel(vbsptlinkagedata.vbspt_trackIDs);
for i=1:n_vbsptTracks
    trackID = vbsptlinkagedata.vbspt_trackIDs(i);
    non_interp_vbspt = vbsptdata.Wbest.est2.sMaxP{i};
    observed_frames = expdata.smLinked{trackID}(:,end) - ...
                      expdata.smLinked{trackID}(1,end) + 1;
    % Final frame does not have a diffusive state.
    expdata.vbspt{trackID} = non_interp_vbspt(observed_frames(1:end-1));
    expdata.vbspt_interp{trackID} = non_interp_vbspt;
end

n_usable_tracks = sum(~toReject);
% note that a specific nMinTracks was used to determine all analprops.
n_clustered_tracks = analprops.totalTracksInClust;

overall.frac_tracks_in_clus = n_clustered_tracks/n_usable_tracks;
assert(numel(vbsptdata.Wbest.est.Ptot)==2,'a non-two state HMM was input.')
overall.slowFrac = vbsptdata.Wbest.est.Ptot(1);

usable_tracks = find(~toReject);
% The below assignment is over-complicated because some of the trackIDs were rotated the wrong way in the cell array [1xN and Nx1]
clustered_tracks=cell2mat(cellfun(@(x) x(:),measprops.trackIDs(analprops.clusters_to_use),'uniformoutput',false));
allnonclusteredtracks = setdiff(usable_tracks,clustered_tracks);

all_slow_frames = 0;
all_slow_frames_interp = 0;
clus_slow_frames = 0;
clus_slow_frames_interp = 0;
nonclus_slow_frames = 0;
nonclus_slow_frames_interp = 0;

all_frames = 0;
all_frames_interp = 0;
clus_frames = 0;
clus_frames_interp = 0;
nonclus_frames = 0;
nonclus_frames_interp = 0;

for i=1:numel(expdata.smLinked)
    if any(i==usable_tracks)
        all_states = expdata.vbspt{i};
        all_states_interp = expdata.vbspt_interp{i};

        all_slow_frames        = all_slow_frames        + sum(all_states==1);
        all_slow_frames_interp = all_slow_frames_interp + sum(all_states_interp==1);
        all_frames             = all_frames             + numel(all_states);
        all_frames_interp      = all_frames_interp      + numel(all_states_interp);

        if any(i==clustered_tracks)
            clus_slow_frames        = clus_slow_frames        + sum(all_states==1);
            clus_slow_frames_interp = clus_slow_frames_interp + sum(all_states_interp==1);
            clus_frames             = clus_frames             + numel(all_states);
            clus_frames_interp      = clus_frames_interp      + numel(all_states_interp);
        end
        if any(i==allnonclusteredtracks)
            nonclus_slow_frames        = nonclus_slow_frames        + sum(all_states==1);
            nonclus_slow_frames_interp = nonclus_slow_frames_interp + sum(all_states_interp==1);
            nonclus_frames             = nonclus_frames             + numel(all_states);
            nonclus_frames_interp      = nonclus_frames_interp      + numel(all_states_interp);
        end
    end
end
overall.nonclus_slowFrac = nonclus_slow_frames/nonclus_frames;
overall.clus_slowFrac = clus_slow_frames/clus_frames;
overall.all_slowFrac = all_slow_frames/all_frames;
overall.nonclus_slowFrac_interp = nonclus_slow_frames_interp/nonclus_frames_interp;
overall.clus_slowFrac_interp = clus_slow_frames_interp/clus_frames_interp;
overall.all_slowFrac_interp = all_slow_frames_interp/all_frames_interp;

% Determining tracks starting in clusters. We repeat this for each definition of cluster region we have.
%% THESE ARE CALCULATED USING THE MASK BOUNDARY DEFINITION!
nclus = numel(measprops.N);
KDEr50_fracRecToClus_maskBound = nan(nclus,1);
KDEr95_fracRecToClus_maskBound = nan(nclus,1);
KDEr50_nRecToClus_maskBound = nan(nclus,1);
KDEr95_nRecToClus_maskBound = nan(nclus,1);
KDEr50_nInClus_maskBound = nan(nclus,1);
KDEr95_nInClus_maskBound = nan(nclus,1);
for i=1:nclus
    % if this cluster was rejected at cluster analysis, we won't work with it here.
    if ~analprops.clusters_to_use(i)
        continue;
    else %this cluster is OK to use
        trackIDs = measprops.trackIDs{i};

        % coordinates of first frame for each track, in um
        first_coords = cell2mat(cellfun(@(x) x(1,1:2)*0.1067,expdata.smLinked(trackIDs),'uniformoutput',false));

        %% kernel density estimate of all localizations, 50% of max density based
		    warning('off')
        cluster_region = polyshape(measprops.allLoc.KDEbound50{i});
        warning('on')

        % check if they are within the cluster.
		    is_rec_to_clus = isinterior(cluster_region,first_coords);
        KDEr50_nRecToClus_maskBound(i) = sum(is_rec_to_clus);
        KDEr50_nInClus_maskBound(i) = numel(trackIDs);
        KDEr50_fracRecToClus_maskBound(i) = KDEr50_nRecToClus_maskBound(i)/KDEr50_nInClus_maskBound(i);

        %% kernel density estimate of all localizations, 95% of max density based
        warning('off')
		    cluster_region = polyshape(measprops.allLoc.KDEbound95{i});
        warning('on')

        % check if they are within the cluster.
		    is_rec_to_clus = isinterior(cluster_region,first_coords);
        KDEr95_nRecToClus_maskBound(i) = sum(is_rec_to_clus);
        KDEr95_nInClus_maskBound(i) = numel(trackIDs);
        KDEr95_fracRecToClus_maskBound(i) = KDEr95_nRecToClus_maskBound(i)/KDEr95_nInClus_maskBound(i);
    end
end
% drop out rejected clusters
KDEr50_toDrop_maskBound  = isnan(KDEr50_fracRecToClus_maskBound);
KDEr95_toDrop_maskBound  = isnan(KDEr95_fracRecToClus_maskBound);

KDEr50_fracRecToClus_maskBound = KDEr50_fracRecToClus_maskBound(~KDEr50_toDrop_maskBound);
KDEr95_fracRecToClus_maskBound = KDEr95_fracRecToClus_maskBound(~KDEr95_toDrop_maskBound);

KDEr50_nRecToClus_maskBound = KDEr50_nRecToClus_maskBound(~KDEr50_toDrop_maskBound);
KDEr95_nRecToClus_maskBound = KDEr95_nRecToClus_maskBound(~KDEr95_toDrop_maskBound);

KDEr50_nInClus_maskBound = KDEr50_nInClus_maskBound(~KDEr50_toDrop_maskBound);
KDEr95_nInClus_maskBound = KDEr95_nInClus_maskBound(~KDEr95_toDrop_maskBound);

overall.KDEr50_fracRecToClus_maskBound.distr = KDEr50_fracRecToClus_maskBound;
overall.KDEr50_fracRecToClus_maskBound.mean = mean(KDEr50_fracRecToClus_maskBound);
overall.KDEr95_fracRecToClus_maskBound.distr = KDEr95_fracRecToClus_maskBound;
overall.KDEr95_fracRecToClus_maskBound.mean = mean(KDEr95_fracRecToClus_maskBound);

overall.KDEr50_nRecToClus_maskBound.distr = KDEr50_nRecToClus_maskBound;
overall.KDEr50_nRecToClus_maskBound.mean = mean(KDEr50_nRecToClus_maskBound);
overall.KDEr95_nRecToClus_maskBound.distr = KDEr95_nRecToClus_maskBound;
overall.KDEr95_nRecToClus_maskBound.mean = mean(KDEr95_nRecToClus_maskBound);

overall.KDEr50_nInClus_maskBound.distr = KDEr50_nInClus_maskBound;
overall.KDEr50_nInClus_maskBound.mean = mean(KDEr50_nInClus_maskBound);
overall.KDEr95_nInClus_maskBound.distr = KDEr95_nInClus_maskBound;
overall.KDEr95_nInClus_maskBound.mean = mean(KDEr95_nInClus_maskBound);

overall.KDEr50_fracRecToClusPooled_maskBound = sum(KDEr50_nRecToClus_maskBound)/sum(KDEr50_nInClus_maskBound);
overall.KDEr95_fracRecToClusPooled_maskBound = sum(KDEr95_nRecToClus_maskBound)/sum(KDEr95_nInClus_maskBound);

% Determine the on/off times of tracks in clusters, and the out-of-cluster distances.
%% THESE ARE CALCULATED USING THE MASK BOUNDARY DEFINITION!

KDEr50_onClus_maskBound = cell(nclus,1);
KDEr95_onClus_maskBound = cell(nclus,1);

KDEr50_outOfClus_dists_maskBound = cell(nclus,1);
KDEr95_outOfClus_dists_maskBound = cell(nclus,1);

KDEr50_outOfClus_startdists_maskBound = cell(nclus,1);
KDEr95_outOfClus_startdists_maskBound = cell(nclus,1);

KDEr50_outOfClus_enddists_maskBound = cell(nclus,1);
KDEr95_outOfClus_enddists_maskBound = cell(nclus,1);

KDEr50_all_dists_maskBound = cell(nclus,1);
KDEr95_all_dists_maskBound = cell(nclus,1);
KDEr50_all_startdists_maskBound = cell(nclus,1);
KDEr95_all_startdists_maskBound = cell(nclus,1);
KDEr50_all_enddists_maskBound = cell(nclus,1);
KDEr95_all_enddists_maskBound = cell(nclus,1);

for i=1:nclus % annotate each track in clusters whether they're in their cluster or not
    % if this cluster was rejected at cluster analysis, we won't work with it here.
    if ~analprops.clusters_to_use(i)
        continue;
    else %this cluster is OK to use
        trackIDs = measprops.trackIDs{i};
        ntracks = numel(trackIDs);

        KDEr50_onClus_maskBound{i} = cell(ntracks,1);
        KDEr95_onClus_maskBound{i} = cell(ntracks,1);

        KDEr50_outOfClus_dists_maskBound{i} = cell(ntracks,1);
        KDEr95_outOfClus_dists_maskBound{i} = cell(ntracks,1);

        KDEr50_outOfClus_startdists_maskBound{i} = cell(ntracks,1);
        KDEr95_outOfClus_startdists_maskBound{i} = cell(ntracks,1);

        KDEr50_outOfClus_enddists_maskBound{i} = cell(ntracks,1);
        KDEr95_outOfClus_enddists_maskBound{i} = cell(ntracks,1);

        KDEr50_all_dists_maskBound{i} = cell(ntracks,1);
        KDEr95_all_dists_maskBound{i} = cell(ntracks,1);
        KDEr50_all_startdists_maskBound{i} = cell(ntracks,1);
        KDEr95_all_startdists_maskBound{i} = cell(ntracks,1);
        KDEr50_all_enddists_maskBound{i} = cell(ntracks,1);
        KDEr95_all_enddists_maskBound{i} = cell(ntracks,1);

        for j=1:ntracks % for each track
            % coordinates for each track, in um
            track_coords = expdata.smLinked{trackIDs(j)}(:,1:2)*0.1067;

            %% kernel density estimate of all localizations, 50% of max density based
            cluster_center = measprops.allLoc.KDEmaskCOM50(i,:);
            warning('off')
			      cluster_region = polyshape(measprops.allLoc.KDEbound50{i});
            warning('off')
            dsq = (track_coords(:,1)-cluster_center(1)).^2 + (track_coords(:,2)-cluster_center(2)).^2;
            KDEr50_all_dists_maskBound{i}{j} = sqrt(dsq);
            KDEr50_all_startdists_maskBound{i}{j} = sqrt(dsq(1));
            KDEr50_all_enddists_maskBound{i}{j} = sqrt(dsq(end));

			      KDEr50_onClus_maskBound{i}{j} = isinterior(cluster_region,track_coords);
            KDEr50_outOfClus_dists_maskBound{i}{j} = sqrt(dsq(~KDEr50_onClus_maskBound{i}{j}));
            if any(~KDEr50_onClus_maskBound{i}{j})
                KDEr50_outOfClus_startdists_maskBound{i}{j} = sqrt(dsq(1));
                KDEr50_outOfClus_enddists_maskBound{i}{j}   = sqrt(dsq(end));
            end

            %% kernel density estimate of all localizations, 95% of max density based
            cluster_center = measprops.allLoc.KDEmaskCOM95(i,:);
            warning('off')
			      cluster_region = polyshape(measprops.allLoc.KDEbound95{i});
            warning('on')
            dsq = (track_coords(:,1)-cluster_center(1)).^2 + (track_coords(:,2)-cluster_center(2)).^2;
            KDEr95_all_dists_maskBound{i}{j} = sqrt(dsq);
            KDEr95_all_startdists_maskBound{i}{j} = sqrt(dsq(1));
            KDEr95_all_enddists_maskBound{i}{j} = sqrt(dsq(end));

            KDEr95_onClus_maskBound{i}{j} = isinterior(cluster_region,track_coords);
            KDEr95_outOfClus_dists_maskBound{i}{j} = sqrt(dsq(~KDEr95_onClus_maskBound{i}{j}));
            if any(~KDEr95_onClus_maskBound{i}{j})
                KDEr95_outOfClus_startdists_maskBound{i}{j} = sqrt(dsq(1));
                KDEr95_outOfClus_enddists_maskBound{i}{j}   = sqrt(dsq(end));
            end
        end
    end
end
overall.KDEr50_outOfClus_dists_maskBound = KDEr50_outOfClus_dists_maskBound(cellfun(@(x) ~isempty(x),KDEr50_outOfClus_dists_maskBound));
overall.KDEr95_outOfClus_dists_maskBound = KDEr95_outOfClus_dists_maskBound(cellfun(@(x) ~isempty(x),KDEr95_outOfClus_dists_maskBound));

overall.KDEr50_outOfClus_startdists_maskBound = KDEr50_outOfClus_startdists_maskBound(cellfun(@(x) ~isempty(x),KDEr50_outOfClus_startdists_maskBound));
overall.KDEr50_outOfClus_enddists_maskBound = KDEr50_outOfClus_enddists_maskBound(cellfun(@(x) ~isempty(x),KDEr50_outOfClus_enddists_maskBound));
overall.KDEr95_outOfClus_startdists_maskBound = KDEr95_outOfClus_startdists_maskBound(cellfun(@(x) ~isempty(x),KDEr95_outOfClus_startdists_maskBound));
overall.KDEr95_outOfClus_enddists_maskBound = KDEr95_outOfClus_enddists_maskBound(cellfun(@(x) ~isempty(x),KDEr95_outOfClus_enddists_maskBound));

overall.KDEr50_all_dists_maskBound = KDEr50_all_dists_maskBound(cellfun(@(x) ~isempty(x),KDEr50_all_dists_maskBound));
overall.KDEr95_all_dists_maskBound = KDEr95_all_dists_maskBound(cellfun(@(x) ~isempty(x),KDEr95_all_dists_maskBound));
overall.KDEr50_all_startdists_maskBound = KDEr50_all_startdists_maskBound(cellfun(@(x) ~isempty(x),KDEr50_all_startdists_maskBound));
overall.KDEr95_all_startdists_maskBound = KDEr95_all_startdists_maskBound(cellfun(@(x) ~isempty(x),KDEr95_all_startdists_maskBound));
overall.KDEr50_all_enddists_maskBound = KDEr50_all_enddists_maskBound(cellfun(@(x) ~isempty(x),KDEr50_all_enddists_maskBound));
overall.KDEr95_all_enddists_maskBound = KDEr95_all_enddists_maskBound(cellfun(@(x) ~isempty(x),KDEr95_all_enddists_maskBound));

% For each cluster provide a cell array of per-track total on time, per-track total off time,
% and then a cell array of cell arrays for per-track continuous on-times and per-track continuous off-times
overall.KDEr50_onClusTime_maskBound.perTrackTotalOnTime = cell(nclus,1);
overall.KDEr95_onClusTime_maskBound.perTrackTotalOnTime = cell(nclus,1);

overall.KDEr50_onClusTime_maskBound.perTrackTotalOffTime = cell(nclus,1);
overall.KDEr95_onClusTime_maskBound.perTrackTotalOffTime = cell(nclus,1);

overall.KDEr50_onClusTime_maskBound.perTrackContinuousOnTime = cell(nclus,1);
overall.KDEr95_onClusTime_maskBound.perTrackContinuousOnTime = cell(nclus,1);

overall.KDEr50_onClusTime_maskBound.perTrackContinuousOffTime = cell(nclus,1);
overall.KDEr95_onClusTime_maskBound.perTrackContinuousOffTime = cell(nclus,1);

for i=1:nclus
    if ~analprops.clusters_to_use(i)
        continue;
    else %this cluster is OK to use
        trackIDs = measprops.trackIDs{i};
        ntracks = numel(trackIDs);

		    overall.KDEr50_onClusTime_maskBound.perTrackTotalOnTime{i} = zeros(ntracks,1);
        overall.KDEr95_onClusTime_maskBound.perTrackTotalOnTime{i} = zeros(ntracks,1);

		    overall.KDEr50_onClusTime_maskBound.perTrackTotalOffTime{i} = zeros(ntracks,1);
        overall.KDEr95_onClusTime_maskBound.perTrackTotalOffTime{i} = zeros(ntracks,1);

		    overall.KDEr50_onClusTime_maskBound.perTrackContinuousOnTime{i} = cell(ntracks,1);
        overall.KDEr95_onClusTime_maskBound.perTrackContinuousOnTime{i} = cell(ntracks,1);

		    overall.KDEr50_onClusTime_maskBound.perTrackContinuousOffTime{i} = cell(ntracks,1);
        overall.KDEr95_onClusTime_maskBound.perTrackContinuousOffTime{i} = cell(ntracks,1);
        for j=1:ntracks
            % the total on cluster, off cluster times for each track are scalars.
            overall.KDEr50_onClusTime_maskBound.perTrackTotalOnTime{i}(j)   = sum( KDEr50_onClus_maskBound{i}{j})*0.02;
            overall.KDEr50_onClusTime_maskBound.perTrackTotalOffTime{i}(j)  = sum(~KDEr50_onClus_maskBound{i}{j})*0.02;
            overall.KDEr95_onClusTime_maskBound.perTrackTotalOnTime{i}(j)   = sum( KDEr95_onClus_maskBound{i}{j})*0.02;
            overall.KDEr95_onClusTime_maskBound.perTrackTotalOffTime{i}(j)  = sum(~KDEr95_onClus_maskBound{i}{j})*0.02;

            % the continuous on cluster, off cluster times for each track are vectors.
            overall.KDEr50_onClusTime_maskBound.perTrackContinuousOnTime{i}{j}   = get_durations_of_true( KDEr50_onClus_maskBound{i}{j}')*0.02;
            overall.KDEr50_onClusTime_maskBound.perTrackContinuousOffTime{i}{j}  = get_durations_of_true(~KDEr50_onClus_maskBound{i}{j}')*0.02;
            overall.KDEr95_onClusTime_maskBound.perTrackContinuousOnTime{i}{j}   = get_durations_of_true( KDEr95_onClus_maskBound{i}{j}')*0.02;
            overall.KDEr95_onClusTime_maskBound.perTrackContinuousOffTime{i}{j}  = get_durations_of_true(~KDEr95_onClus_maskBound{i}{j}')*0.02;
        end

    end
end

end

function durations = get_durations_of_true(boolvec)
if ~islogical(boolvec)
    error('input a logical vector')
end
boolvec=~boolvec; % the following lines of code actually find durations of ZEROS.
dbvec = diff([true boolvec true]);
startIDX = find(dbvec<0);
endIDX   = find(dbvec>0)-1;
durations = endIDX-startIDX+1;
end

function [all_path,all_file] = get_file_structures(datafolders)
%% Get file structures (lvPALM files)
%datafolders = {'/Volumes/med/Hahn_Lab/Lab Members/Bei/Data Backup/Binder-tag/03282018 - SLBH/100pg/'};
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
end

function [vbSPTfile,vbSPTmetafile] = get_vbSPT_file_structures(ncells,vbSPTdata_dir,vbSPTmeta_dir)
vbSPTfile.tag = cell(ncells,1);
vbSPTfile.bin = cell(ncells,1);
vbSPTmetafile.tag = cell(ncells,1);
vbSPTmetafile.bin = cell(ncells,1);

%vbSPTdata_dir = './vbspt/results/';
%vbSPTmeta_dir = './vbspt/vbspt_source/';
for i=1:ncells
    vbSPTfile.tag{i} = sprintf('%stag_%02d_HMManalysis_hidden2.mat',vbSPTdata_dir,i);
    vbSPTfile.bin{i} = sprintf('%sbin_%02d_HMManalysis_hidden2.mat',vbSPTdata_dir,i);
    vbSPTmetafile.tag{i} = sprintf('%stag_%02d.mat',vbSPTmeta_dir,i);
    vbSPTmetafile.bin{i} = sprintf('%sbin_%02d.mat',vbSPTmeta_dir,i);
end


end
