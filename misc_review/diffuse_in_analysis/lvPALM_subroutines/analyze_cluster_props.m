% ANALYZE_CLUSTER_PROPS
% Analyzes and visualizes cluster properties, including various estimates of cluster size,
% cluster region, # tracks per cluster, slow and fast diffusive zones...
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
%   SPTacqTime: duration of SPT data acquisition block (seconds)
%
%   RIcutoff: cutoff used for filtering recruitment intervals
%
% Saves a .mat file to the working directory:
%   [dateprefix filesuffix '.mat']
%
% Also produces figures for each type of analysis carried out.
%
% Part of the cluster_segmentation.m pipeline.
% 2018 May 9 / Mike Pablo

function analyze_cluster_props(dateprefix,propsuffix,filesuffix,nMinTracks,RIcutoff,parentFigDir)
clusprop = load([dateprefix propsuffix]);
ncells = numel(clusprop.tag_prop);
% propPointEst: point estimates of the property distribution
% Can be the arithmetic mean, lognormal fitted mean, exponential fitted half-life, etc.
% See the analyze_props function to determine how each property point estimate was obtained.
propPointEst_tag = cell(ncells,1);
propPointEst_bin = cell(ncells,1);
propdists_tag = cell(ncells,1);
propdists_bin = cell(ncells,1);

figSubDir = 'clus_analysis/';
parentFigDir = [parentFigDir figSubDir];
warning('off');
mkdir(parentFigDir);
warning('on');

parfor i=1:ncells
    [propPointEst_tag{i},propdists_tag{i}] = analyze_props(clusprop.tag_prop{i},nMinTracks);
    [propPointEst_bin{i},propdists_bin{i}] = analyze_props(clusprop.bin_prop{i},nMinTracks);
end
% for compatibility with pipeline (may 11 2018) save results in different name
propPointEst.tag = propPointEst_tag;
propdists.tag = propdists_tag;
propPointEst.bin = propPointEst_bin;
propdists.bin = propdists_bin;
save([dateprefix filesuffix],'propPointEst','propdists','nMinTracks','RIcutoff');
end


function [propPointEst,propdistr] = analyze_props(data,nMinTracks)
disp('Calculating mean properties and property distributions.')

clusters_to_use = (data.N >= nMinTracks);

% store the clusters_to_use parameter under propPoint convenience.
propPointEst.clusters_to_use = clusters_to_use;

% track numbers, arithmetic mean
propPointEst.N = mean(data.N(clusters_to_use));
propdistr.N = data.N(clusters_to_use);

propPointEst.totalTracksInClust = sum(data.N(clusters_to_use));

% whether multi-region segmentation was observed in the KDE
propdistr.allLoc.KDE.multiRegion50 = ~isnan(data.allLoc.KDEeffR50(clusters_to_use)) & isnan(data.allLoc.KDEcirc50(clusters_to_use));
propPointEst.allLoc.KDE.multiRegion50 = sum(propdistr.allLoc.KDE.multiRegion50)/numel(propdistr.allLoc.KDE.multiRegion50);
propdistr.allLoc.KDE.multiRegion95 = ~isnan(data.allLoc.KDEeffR95(clusters_to_use)) & isnan(data.allLoc.KDEcirc95(clusters_to_use));
propPointEst.allLoc.KDE.multiRegion95 = sum(propdistr.allLoc.KDE.multiRegion95)/numel(propdistr.allLoc.KDE.multiRegion95);

% % cluster circularity
propdistr.allLoc.KDE.circ50 = data.allLoc.KDEcirc50(clusters_to_use);
propPointEst.allLoc.KDE.circ50 = nanmedian(propdistr.allLoc.KDE.circ50); % value is nan if KDE attempted with 0 or 1 points,
                                                                       % OR if multiple regions were found at thresholding
propdistr.allLoc.KDE.circ95 = data.allLoc.KDEcirc95(clusters_to_use);
propPointEst.allLoc.KDE.circ95 = nanmedian(propdistr.allLoc.KDE.circ95);

% cluster eccentricity
propdistr.allLoc.KDE.eccen50 = data.allLoc.KDEeccen50(clusters_to_use);
propPointEst.allLoc.KDE.eccen50 = nanmedian(propdistr.allLoc.KDE.eccen50); % value is nan if KDE attempted with 0 or 1 points,
                                                                       % OR if multiple regions were found at thresholding
propdistr.allLoc.KDE.eccen95 = data.allLoc.KDEeccen95(clusters_to_use);
propPointEst.allLoc.KDE.eccen95 = nanmedian(propdistr.allLoc.KDE.eccen95);

% cluster area
propdistr.allLoc.KDE.area50 = data.allLoc.KDEarea50(clusters_to_use);
propPointEst.allLoc.KDE.area50 = nanmedian(propdistr.allLoc.KDE.area50);

propdistr.allLoc.KDE.area95 = data.allLoc.KDEarea95(clusters_to_use);
propPointEst.allLoc.KDE.area95 = nanmedian(propdistr.allLoc.KDE.area95);

% cluster major axis
propdistr.allLoc.KDE.majAx50 = data.allLoc.KDEmajAx50(clusters_to_use);
propPointEst.allLoc.KDE.majAx50 = nanmedian(propdistr.allLoc.KDE.majAx50);

propdistr.allLoc.KDE.majAx95 = data.allLoc.KDEmajAx95(clusters_to_use);
propPointEst.allLoc.KDE.majAx95 = nanmedian(propdistr.allLoc.KDE.majAx95);

% cluster minor axis
propdistr.allLoc.KDE.minAx50 = data.allLoc.KDEminAx50(clusters_to_use);
propPointEst.allLoc.KDE.minAx50 = nanmedian(propdistr.allLoc.KDE.minAx50);

propdistr.allLoc.KDE.minAx95 = data.allLoc.KDEminAx95(clusters_to_use);
propPointEst.allLoc.KDE.minAx95 = nanmedian(propdistr.allLoc.KDE.minAx95);

% cluster boundaries
propdistr.allLoc.KDE.bound50 = data.allLoc.KDEbound50(clusters_to_use);
propdistr.allLoc.KDE.bound95 = data.allLoc.KDEbound95(clusters_to_use);

% cluster mask COMs
propdistr.allLoc.KDEmaskCOM50 = data.allLoc.KDEmaskCOM50(clusters_to_use,:);
propdistr.allLoc.KDEmaskCOM95 = data.allLoc.KDEmaskCOM95(clusters_to_use,:);

% centroid-based cluster size, median
propdistr.cent.fitR = data.cent.fitR(clusters_to_use);
propPointEst.cent.fitR = nanmedian(propdistr.cent.fitR);

% slow fraction, arithmetic mean
propPointEst.vbSPT.slowFrac = mean(data.vbSPT.slowFrac(clusters_to_use));
propdistr.vbSPT.slowFrac = data.vbSPT.slowFrac(clusters_to_use);

% % other fast and slow properties from KDE analysis
propdistr.vbSPT.slowcirc50 = data.vbSPT.slowcirc50(clusters_to_use);
propPointEst.vbSPT.slowcirc50 = nanmedian(propdistr.vbSPT.slowcirc50);
propdistr.vbSPT.slowcirc95 = data.vbSPT.slowcirc95(clusters_to_use);
propPointEst.vbSPT.slowcirc95 = nanmedian(propdistr.vbSPT.slowcirc95);
propdistr.vbSPT.fastcirc50 = data.vbSPT.fastcirc50(clusters_to_use);
propPointEst.vbSPT.fastcirc50 = nanmedian(propdistr.vbSPT.fastcirc50);
propdistr.vbSPT.fastcirc95 = data.vbSPT.fastcirc95(clusters_to_use);
propPointEst.vbSPT.fastcirc95 = nanmedian(propdistr.vbSPT.fastcirc95);

propdistr.vbSPT.sloweccen50 = data.vbSPT.sloweccen50(clusters_to_use);
propPointEst.vbSPT.sloweccen50 = nanmedian(propdistr.vbSPT.sloweccen50);
propdistr.vbSPT.sloweccen95 = data.vbSPT.sloweccen95(clusters_to_use);
propPointEst.vbSPT.sloweccen95 = nanmedian(propdistr.vbSPT.sloweccen95);
propdistr.vbSPT.fasteccen50 = data.vbSPT.fasteccen50(clusters_to_use);
propPointEst.vbSPT.fasteccen50 = nanmedian(propdistr.vbSPT.fasteccen50);
propdistr.vbSPT.fasteccen95 = data.vbSPT.fasteccen95(clusters_to_use);
propPointEst.vbSPT.fasteccen95 = nanmedian(propdistr.vbSPT.fasteccen95);

propdistr.vbSPT.slowarea50 = data.vbSPT.slowarea50(clusters_to_use);
propPointEst.vbSPT.slowarea50 = nanmedian(propdistr.vbSPT.slowarea50);
propdistr.vbSPT.slowarea95 = data.vbSPT.slowarea95(clusters_to_use);
propPointEst.vbSPT.slowarea95 = nanmedian(propdistr.vbSPT.slowarea95);
propdistr.vbSPT.fastarea50 = data.vbSPT.fastarea50(clusters_to_use);
propPointEst.vbSPT.fastarea50 = nanmedian(propdistr.vbSPT.fastarea50);
propdistr.vbSPT.fastarea95 = data.vbSPT.fastarea95(clusters_to_use);
propPointEst.vbSPT.fastarea95 = nanmedian(propdistr.vbSPT.fastarea95);

propdistr.vbSPT.slowmajAx50 = data.vbSPT.slowmajAx50(clusters_to_use);
propPointEst.vbSPT.slowmajAx50 = nanmedian(propdistr.vbSPT.slowmajAx50);
propdistr.vbSPT.slowmajAx95 = data.vbSPT.slowmajAx95(clusters_to_use);
propPointEst.vbSPT.slowmajAx95 = nanmedian(propdistr.vbSPT.slowmajAx95);
propdistr.vbSPT.fastmajAx50 = data.vbSPT.fastmajAx50(clusters_to_use);
propPointEst.vbSPT.fastmajAx50 = nanmedian(propdistr.vbSPT.fastmajAx50);
propdistr.vbSPT.fastmajAx95 = data.vbSPT.fastmajAx95(clusters_to_use);
propPointEst.vbSPT.fastmajAx95 = nanmedian(propdistr.vbSPT.fastmajAx95);

propdistr.vbSPT.slowminAx50 = data.vbSPT.slowminAx50(clusters_to_use);
propPointEst.vbSPT.slowminAx50 = nanmedian(propdistr.vbSPT.slowminAx50);
propdistr.vbSPT.slowminAx95 = data.vbSPT.slowminAx95(clusters_to_use);
propPointEst.vbSPT.slowminAx95 = nanmedian(propdistr.vbSPT.slowminAx95);
propdistr.vbSPT.fastminAx50 = data.vbSPT.fastminAx50(clusters_to_use);
propPointEst.vbSPT.fastminAx50 = nanmedian(propdistr.vbSPT.fastminAx50);
propdistr.vbSPT.fastminAx95 = data.vbSPT.fastminAx95(clusters_to_use);
propPointEst.vbSPT.fastminAx95 = nanmedian(propdistr.vbSPT.fastminAx95);

propdistr.vbSPT.slowbound50 = data.vbSPT.slowbound50(clusters_to_use);
propdistr.vbSPT.slowbound95 = data.vbSPT.slowbound95(clusters_to_use);
propdistr.vbSPT.fastbound50 = data.vbSPT.fastbound50(clusters_to_use);
propdistr.vbSPT.fastbound95 = data.vbSPT.fastbound95(clusters_to_use);

propdistr.vbSPT.slowCOM50 = data.vbSPT.slowCOM50(clusters_to_use,:);
propdistr.vbSPT.slowCOM95 = data.vbSPT.slowCOM95(clusters_to_use,:);
propdistr.vbSPT.fastCOM50 = data.vbSPT.fastCOM50(clusters_to_use,:);
propdistr.vbSPT.fastCOM95 = data.vbSPT.fastCOM95(clusters_to_use,:);

% [REMOVE DIFFUSE IN RECALCULATIONS]===================================
propdistr.allLoc.KDE.multiRegion50_removeDiffuseIn = ~isnan(data.allLoc.KDEeffR50_removeDiffuseIn(clusters_to_use)) & isnan(data.allLoc.KDEcirc50_removeDiffuseIn(clusters_to_use));
propPointEst.allLoc.KDE.multiRegion50_removeDiffuseIn = sum(propdistr.allLoc.KDE.multiRegion50_removeDiffuseIn)/numel(propdistr.allLoc.KDE.multiRegion50_removeDiffuseIn);
propdistr.allLoc.KDE.multiRegion95_removeDiffuseIn = ~isnan(data.allLoc.KDEeffR95_removeDiffuseIn(clusters_to_use)) & isnan(data.allLoc.KDEcirc95_removeDiffuseIn(clusters_to_use));
propPointEst.allLoc.KDE.multiRegion95_removeDiffuseIn = sum(propdistr.allLoc.KDE.multiRegion95_removeDiffuseIn)/numel(propdistr.allLoc.KDE.multiRegion95_removeDiffuseIn);

% % cluster circularity
propdistr.allLoc.KDE.circ50_removeDiffuseIn = data.allLoc.KDEcirc50_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.circ50_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.circ50_removeDiffuseIn); % value is nan if KDE attempted with 0 or 1 points,
                                                                       % OR if multiple regions were found at thresholding
propdistr.allLoc.KDE.circ95_removeDiffuseIn = data.allLoc.KDEcirc95_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.circ95_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.circ95_removeDiffuseIn);

% cluster eccentricity
propdistr.allLoc.KDE.eccen50_removeDiffuseIn = data.allLoc.KDEeccen50_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.eccen50_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.eccen50_removeDiffuseIn); % value is nan if KDE attempted with 0 or 1 points,
                                                                       % OR if multiple regions were found at thresholding
propdistr.allLoc.KDE.eccen95_removeDiffuseIn = data.allLoc.KDEeccen95_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.eccen95_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.eccen95_removeDiffuseIn);

% cluster area
propdistr.allLoc.KDE.area50_removeDiffuseIn = data.allLoc.KDEarea50_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.area50_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.area50_removeDiffuseIn);

propdistr.allLoc.KDE.area95_removeDiffuseIn = data.allLoc.KDEarea95_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.area95_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.area95_removeDiffuseIn);

% cluster major axis
propdistr.allLoc.KDE.majAx50_removeDiffuseIn = data.allLoc.KDEmajAx50_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.majAx50_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.majAx50_removeDiffuseIn);

propdistr.allLoc.KDE.majAx95_removeDiffuseIn = data.allLoc.KDEmajAx95_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.majAx95_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.majAx95_removeDiffuseIn);

% cluster minor axis
propdistr.allLoc.KDE.minAx50_removeDiffuseIn = data.allLoc.KDEminAx50_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.minAx50_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.minAx50_removeDiffuseIn);

propdistr.allLoc.KDE.minAx95_removeDiffuseIn = data.allLoc.KDEminAx95_removeDiffuseIn(clusters_to_use);
propPointEst.allLoc.KDE.minAx95_removeDiffuseIn = nanmedian(propdistr.allLoc.KDE.minAx95_removeDiffuseIn);

% cluster boundaries
propdistr.allLoc.KDE.bound50_removeDiffuseIn = data.allLoc.KDEbound50_removeDiffuseIn(clusters_to_use);
propdistr.allLoc.KDE.bound95_removeDiffuseIn = data.allLoc.KDEbound95_removeDiffuseIn(clusters_to_use);

% cluster mask COMs
propdistr.allLoc.KDEmaskCOM50_removeDiffuseIn = data.allLoc.KDEmaskCOM50_removeDiffuseIn(clusters_to_use,:);
propdistr.allLoc.KDEmaskCOM95_removeDiffuseIn = data.allLoc.KDEmaskCOM95_removeDiffuseIn(clusters_to_use,:);


% % other fast and slow properties from KDE analysis
propdistr.vbSPT.slowcirc50_removeDiffuseIn = data.vbSPT.slowcirc50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.slowcirc50_removeDiffuseIn = nanmedian(propdistr.vbSPT.slowcirc50_removeDiffuseIn);
propdistr.vbSPT.slowcirc95_removeDiffuseIn = data.vbSPT.slowcirc95_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.slowcirc95_removeDiffuseIn = nanmedian(propdistr.vbSPT.slowcirc95_removeDiffuseIn);
propdistr.vbSPT.fastcirc50_removeDiffuseIn = data.vbSPT.fastcirc50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fastcirc50_removeDiffuseIn = nanmedian(propdistr.vbSPT.fastcirc50_removeDiffuseIn);
propdistr.vbSPT.fastcirc95_removeDiffuseIn = data.vbSPT.fastcirc95_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fastcirc95_removeDiffuseIn = nanmedian(propdistr.vbSPT.fastcirc9_removeDiffuseIn5);

propdistr.vbSPT.sloweccen50_removeDiffuseIn = data.vbSPT.sloweccen50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.sloweccen50_removeDiffuseIn = nanmedian(propdistr.vbSPT.sloweccen50_removeDiffuseIn);
propdistr.vbSPT.sloweccen9_removeDiffuseIn5 = data.vbSPT.sloweccen95_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.sloweccen95_removeDiffuseIn = nanmedian(propdistr.vbSPT.sloweccen95_removeDiffuseIn);
propdistr.vbSPT.fasteccen50_removeDiffuseIn = data.vbSPT.fasteccen50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fasteccen50_removeDiffuseIn = nanmedian(propdistr.vbSPT.fasteccen50_removeDiffuseIn);
propdistr.vbSPT.fasteccen95_removeDiffuseIn = data.vbSPT.fasteccen95_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fasteccen95_removeDiffuseIn = nanmedian(propdistr.vbSPT.fasteccen95_removeDiffuseIn);

propdistr.vbSPT.slowarea50_removeDiffuseIn = data.vbSPT.slowarea50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.slowarea50_removeDiffuseIn = nanmedian(propdistr.vbSPT.slowarea50_removeDiffuseIn);
propdistr.vbSPT.slowarea95_removeDiffuseIn = data.vbSPT.slowarea9_removeDiffuseIn5(clusters_to_use);
propPointEst.vbSPT.slowarea95_removeDiffuseIn = nanmedian(propdistr.vbSPT.slowarea95_removeDiffuseIn);
propdistr.vbSPT.fastarea50_removeDiffuseIn = data.vbSPT.fastarea50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fastarea50_removeDiffuseIn = nanmedian(propdistr.vbSPT.fastarea50_removeDiffuseIn);
propdistr.vbSPT.fastarea95_removeDiffuseIn = data.vbSPT.fastarea95_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fastarea95_removeDiffuseIn = nanmedian(propdistr.vbSPT.fastarea95_removeDiffuseIn);

propdistr.vbSPT.slowmajAx50_removeDiffuseIn = data.vbSPT.slowmajAx50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.slowmajAx50_removeDiffuseIn = nanmedian(propdistr.vbSPT.slowmajAx50_removeDiffuseIn);
propdistr.vbSPT.slowmajAx95_removeDiffuseIn = data.vbSPT.slowmajAx95_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.slowmajAx95_removeDiffuseIn = nanmedian(propdistr.vbSPT.slowmajAx95_removeDiffuseIn);
propdistr.vbSPT.fastmajAx50_removeDiffuseIn = data.vbSPT.fastmajAx50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fastmajAx50_removeDiffuseIn = nanmedian(propdistr.vbSPT.fastmajAx50_removeDiffuseIn);
propdistr.vbSPT.fastmajAx95_removeDiffuseIn = data.vbSPT.fastmajAx95_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fastmajAx95_removeDiffuseIn = nanmedian(propdistr.vbSPT.fastmajAx95_removeDiffuseIn);

propdistr.vbSPT.slowminAx50_removeDiffuseIn = data.vbSPT.slowminAx50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.slowminAx50_removeDiffuseIn = nanmedian(propdistr.vbSPT.slowminAx50_removeDiffuseIn);
propdistr.vbSPT.slowminAx95_removeDiffuseIn = data.vbSPT.slowminAx95_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.slowminAx95_removeDiffuseIn = nanmedian(propdistr.vbSPT.slowminAx95_removeDiffuseIn);
propdistr.vbSPT.fastminAx50_removeDiffuseIn = data.vbSPT.fastminAx50_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fastminAx50_removeDiffuseIn = nanmedian(propdistr.vbSPT.fastminAx50_removeDiffuseIn);
propdistr.vbSPT.fastminAx95_removeDiffuseIn = data.vbSPT.fastminAx95_removeDiffuseIn(clusters_to_use);
propPointEst.vbSPT.fastminAx95_removeDiffuseIn = nanmedian(propdistr.vbSPT.fastminAx95_removeDiffuseIn);

propdistr.vbSPT.slowbound50_removeDiffuseIn = data.vbSPT.slowbound50_removeDiffuseIn(clusters_to_use);
propdistr.vbSPT.slowbound95_removeDiffuseIn = data.vbSPT.slowbound95_removeDiffuseIn(clusters_to_use);
propdistr.vbSPT.fastbound50_removeDiffuseIn = data.vbSPT.fastbound50_removeDiffuseIn(clusters_to_use);
propdistr.vbSPT.fastbound95_removeDiffuseIn = data.vbSPT.fastbound95_removeDiffuseIn(clusters_to_use);

propdistr.vbSPT.slowCOM50_removeDiffuseIn = data.vbSPT.slowCOM50_removeDiffuseIn(clusters_to_use,:);
propdistr.vbSPT.slowCOM95_removeDiffuseIn = data.vbSPT.slowCOM95_removeDiffuseIn(clusters_to_use,:);
propdistr.vbSPT.fastCOM50_removeDiffuseIn = data.vbSPT.fastCOM50_removeDiffuseIn(clusters_to_use,:);
propdistr.vbSPT.fastCOM95_removeDiffuseIn = data.vbSPT.fastCOM95_removeDiffuseIn(clusters_to_use,:);
% [END REMOVE DIFFUSE IN RECALCULATIONS]===================================

% Lifetime, and exponential fits; the fits were exploratory but left in here.
propdistr.lifetime = data.lifetime(clusters_to_use);

% recruitment intervals
RI_subset = cell2mat(data.rec_int(clusters_to_use));
propdistr.rec_int_by_cluster = data.rec_int(clusters_to_use);
propdistr.rec_int = RI_subset;
end
