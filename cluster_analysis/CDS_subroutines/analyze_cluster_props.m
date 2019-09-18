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
%   RIcutoff: cutoff used for filtering recruitment intervals
%
% Saves a .mat file to the working directory:
%   [dateprefix filesuffix '.mat']
%
% Also produces figures for each type of analysis carried out.
%
% Part of the cluster_segmentation.m pipeline.
% 2018 May 9 / Mike Pablo

function analyze_cluster_props(dateprefix,propsuffix,filesuffix,nMinTracks,RIcutoff)
clusprop = load([dateprefix propsuffix]);
ncells = numel(clusprop.tag_prop);
% propPointEst: point estimates of the property distribution
% Can be the arithmetic mean, lognormal fitted mean, exponential fitted half-life, etc.
% See the analyze_props function to determine how each property point estimate was obtained.
propPointEst_tag = cell(ncells,1);
propPointEst_bin = cell(ncells,1);
propdists_tag = cell(ncells,1);
propdists_bin = cell(ncells,1);
propPointEst_tag_byAd = cell(ncells,1);
propPointEst_bin_byAd = cell(ncells,1);
propdists_tag_byAd = cell(ncells,1);
propdists_bin_byAd = cell(ncells,1);

parfor i=1:ncells
    [propPointEst_tag{i},propdists_tag{i}] = analyze_props(clusprop.tag_prop{i},nMinTracks);
    [propPointEst_bin{i},propdists_bin{i}] = analyze_props(clusprop.bin_prop{i},nMinTracks);

    [propPointEst_tag_byAd{i},propdists_tag_byAd{i}] = analyze_props_by_adhesion(clusprop.tag_prop{i},nMinTracks);
    [propPointEst_bin_byAd{i},propdists_bin_byAd{i}] = analyze_props_by_adhesion(clusprop.bin_prop{i},nMinTracks);
end
% for compatibility with pipeline (may 11 2018) save results in different name
propPointEst.tag = propPointEst_tag;
propdists.tag = propdists_tag;
propPointEst.bin = propPointEst_bin;
propdists.bin = propdists_bin;

propPointEst_byAd.tag = propPointEst_tag_byAd;
propdists_byAd.tag = propdists_tag_byAd;
propPointEst_byAd.bin = propPointEst_bin_byAd;
propdists_byAd.bin = propdists_bin_byAd;

save([dateprefix filesuffix],'propPointEst','propdists','propPointEst_byAd','propdists_byAd',...
                             'nMinTracks','RIcutoff');

end


function [propPointEst,propdistr] = analyze_props_by_adhesion(data,nMinTracks)
disp('Calculating mean properties and property distributions by adhesion association.')

overall_clusters_to_use = (data.N >= nMinTracks);% & (data.cent.fitR<maxClusSize);

% use a 0.1067*2 (2px) cutoff to define whether a cluster is associated with an adhesion.
ad_assoc_thresh = 0.1067*2; % um
ad_associated_clusters = (data.adhesion.d2ad < ad_assoc_thresh) & overall_clusters_to_use;
propPointEst.ad_assoc_type.NA = ad_associated_clusters & (data.adsize'>= 5) & (data.adsize'<22);
propPointEst.ad_assoc_type.FC = ad_associated_clusters & (data.adsize'>=22) & (data.adsize'<88);
propPointEst.ad_assoc_type.FA = ad_associated_clusters & (data.adsize'>=88) & (data.adsize'<440);
propPointEst.ad_assoc_type.off = ~(data.adhesion.d2ad < ad_assoc_thresh) & overall_clusters_to_use;
propPointEst.ad_assoc_type.on  = propPointEst.ad_assoc_type.NA | propPointEst.ad_assoc_type.FC | propPointEst.ad_assoc_type.FA;
propPointEst.ad_assoc_type.small = ad_associated_clusters & (data.adsize'>= 5) & (data.adsize'<88);
propPointEst.ad_assoc_type.large = ad_associated_clusters & (data.adsize'>=88) & (data.adsize'<440);

propPointEst.overall_clusters_to_use = overall_clusters_to_use;

adhesion_labels = {'NA','FC','FA','off','on','small','large'};
% we need to measure all the normal properties, but separated out by adhesion type.

for i=1:numel(adhesion_labels)
    curr_field = adhesion_labels{i};
    clusters_to_use = propPointEst.ad_assoc_type.(curr_field);

    % track numbers, arithmetic mean
    propPointEst.(curr_field).N = mean(data.N(clusters_to_use));
    propdistr.(curr_field).N = data.N(clusters_to_use);

    propPointEst.(curr_field).totalTracksInClust = sum(data.N(clusters_to_use));

    % whether multi-region segmentation was observed in the KDE
    propdistr.(curr_field).allLoc.KDE.multiRegion50 = ~isnan(data.allLoc.KDEeffR50(clusters_to_use)) & isnan(data.allLoc.KDEcirc50(clusters_to_use));
    propPointEst.(curr_field).allLoc.KDE.multiRegion50 = sum(propdistr.(curr_field).allLoc.KDE.multiRegion50)/numel(propdistr.(curr_field).allLoc.KDE.multiRegion50);
    propdistr.(curr_field).allLoc.KDE.multiRegion95 = ~isnan(data.allLoc.KDEeffR95(clusters_to_use)) & isnan(data.allLoc.KDEcirc95(clusters_to_use));
    propPointEst.(curr_field).allLoc.KDE.multiRegion95 = sum(propdistr.(curr_field).allLoc.KDE.multiRegion95)/numel(propdistr.(curr_field).allLoc.KDE.multiRegion95);

    % cluster eccentricity
    propdistr.(curr_field).allLoc.KDE.eccen50 = data.allLoc.KDEeccen50(clusters_to_use);
    propPointEst.(curr_field).allLoc.KDE.eccen50 = nanmedian(propdistr.(curr_field).allLoc.KDE.eccen50); % value is nan if KDE attempted with 0 or 1 points,
                                                                           % OR if multiple regions were found at thresholding
    propdistr.(curr_field).allLoc.KDE.eccen95 = data.allLoc.KDEeccen95(clusters_to_use);
    propPointEst.(curr_field).allLoc.KDE.eccen95 = nanmedian(propdistr.(curr_field).allLoc.KDE.eccen95);

    % cluster area
    propdistr.(curr_field).allLoc.KDE.area50 = data.allLoc.KDEarea50(clusters_to_use);
    propPointEst.(curr_field).allLoc.KDE.area50 = nanmedian(propdistr.(curr_field).allLoc.KDE.area50);

    propdistr.(curr_field).allLoc.KDE.area95 = data.allLoc.KDEarea95(clusters_to_use);
    propPointEst.(curr_field).allLoc.KDE.area95 = nanmedian(propdistr.(curr_field).allLoc.KDE.area95);

    % cluster major axis
    propdistr.(curr_field).allLoc.KDE.majAx50 = data.allLoc.KDEmajAx50(clusters_to_use);
    propPointEst.(curr_field).allLoc.KDE.majAx50 = nanmedian(propdistr.(curr_field).allLoc.KDE.majAx50);

    propdistr.(curr_field).allLoc.KDE.majAx95 = data.allLoc.KDEmajAx95(clusters_to_use);
    propPointEst.(curr_field).allLoc.KDE.majAx95 = nanmedian(propdistr.(curr_field).allLoc.KDE.majAx95);

    % cluster minor axis
    propdistr.(curr_field).allLoc.KDE.minAx50 = data.allLoc.KDEminAx50(clusters_to_use);
    propPointEst.(curr_field).allLoc.KDE.minAx50 = nanmedian(propdistr.(curr_field).allLoc.KDE.minAx50);

    propdistr.(curr_field).allLoc.KDE.minAx95 = data.allLoc.KDEminAx95(clusters_to_use);
    propPointEst.(curr_field).allLoc.KDE.minAx95 = nanmedian(propdistr.(curr_field).allLoc.KDE.minAx95);

    % cluster boundaries
    propdistr.(curr_field).allLoc.KDE.bound50 = data.allLoc.KDEbound50(clusters_to_use);
    propdistr.(curr_field).allLoc.KDE.bound95 = data.allLoc.KDEbound95(clusters_to_use);

    % cluster mask COMs
    propdistr.(curr_field).allLoc.KDEmaskCOM50 = data.allLoc.KDEmaskCOM50(clusters_to_use,:);
    propdistr.(curr_field).allLoc.KDEmaskCOM95 = data.allLoc.KDEmaskCOM95(clusters_to_use,:);

    % centroid-based cluster size, median
    propdistr.(curr_field).cent.fitR = data.cent.fitR(clusters_to_use);
    propPointEst.(curr_field).cent.fitR = nanmedian(propdistr.(curr_field).cent.fitR);

    % slow fraction, arithmetic mean
    propPointEst.(curr_field).vbSPT.slowFrac = mean(data.vbSPT.slowFrac(clusters_to_use));
    propdistr.(curr_field).vbSPT.slowFrac = data.vbSPT.slowFrac(clusters_to_use);

    propdistr.(curr_field).vbSPT.sloweccen50 = data.vbSPT.sloweccen50(clusters_to_use);
    propPointEst.(curr_field).vbSPT.sloweccen50 = nanmedian(propdistr.(curr_field).vbSPT.sloweccen50);
    propdistr.(curr_field).vbSPT.sloweccen95 = data.vbSPT.sloweccen95(clusters_to_use);
    propPointEst.(curr_field).vbSPT.sloweccen95 = nanmedian(propdistr.(curr_field).vbSPT.sloweccen95);
    propdistr.(curr_field).vbSPT.fasteccen50 = data.vbSPT.fasteccen50(clusters_to_use);
    propPointEst.(curr_field).vbSPT.fasteccen50 = nanmedian(propdistr.(curr_field).vbSPT.fasteccen50);
    propdistr.(curr_field).vbSPT.fasteccen95 = data.vbSPT.fasteccen95(clusters_to_use);
    propPointEst.(curr_field).vbSPT.fasteccen95 = nanmedian(propdistr.(curr_field).vbSPT.fasteccen95);

    propdistr.(curr_field).vbSPT.slowarea50 = data.vbSPT.slowarea50(clusters_to_use);
    propPointEst.(curr_field).vbSPT.slowarea50 = nanmedian(propdistr.(curr_field).vbSPT.slowarea50);
    propdistr.(curr_field).vbSPT.slowarea95 = data.vbSPT.slowarea95(clusters_to_use);
    propPointEst.(curr_field).vbSPT.slowarea95 = nanmedian(propdistr.(curr_field).vbSPT.slowarea95);
    propdistr.(curr_field).vbSPT.fastarea50 = data.vbSPT.fastarea50(clusters_to_use);
    propPointEst.(curr_field).vbSPT.fastarea50 = nanmedian(propdistr.(curr_field).vbSPT.fastarea50);
    propdistr.(curr_field).vbSPT.fastarea95 = data.vbSPT.fastarea95(clusters_to_use);
    propPointEst.(curr_field).vbSPT.fastarea95 = nanmedian(propdistr.(curr_field).vbSPT.fastarea95);

    propdistr.(curr_field).vbSPT.slowmajAx50 = data.vbSPT.slowmajAx50(clusters_to_use);
    propPointEst.(curr_field).vbSPT.slowmajAx50 = nanmedian(propdistr.(curr_field).vbSPT.slowmajAx50);
    propdistr.(curr_field).vbSPT.slowmajAx95 = data.vbSPT.slowmajAx95(clusters_to_use);
    propPointEst.(curr_field).vbSPT.slowmajAx95 = nanmedian(propdistr.(curr_field).vbSPT.slowmajAx95);
    propdistr.(curr_field).vbSPT.fastmajAx50 = data.vbSPT.fastmajAx50(clusters_to_use);
    propPointEst.(curr_field).vbSPT.fastmajAx50 = nanmedian(propdistr.(curr_field).vbSPT.fastmajAx50);
    propdistr.(curr_field).vbSPT.fastmajAx95 = data.vbSPT.fastmajAx95(clusters_to_use);
    propPointEst.(curr_field).vbSPT.fastmajAx95 = nanmedian(propdistr.(curr_field).vbSPT.fastmajAx95);

    propdistr.(curr_field).vbSPT.slowminAx50 = data.vbSPT.slowminAx50(clusters_to_use);
    propPointEst.(curr_field).vbSPT.slowminAx50 = nanmedian(propdistr.(curr_field).vbSPT.slowminAx50);
    propdistr.(curr_field).vbSPT.slowminAx95 = data.vbSPT.slowminAx95(clusters_to_use);
    propPointEst.(curr_field).vbSPT.slowminAx95 = nanmedian(propdistr.(curr_field).vbSPT.slowminAx95);
    propdistr.(curr_field).vbSPT.fastminAx50 = data.vbSPT.fastminAx50(clusters_to_use);
    propPointEst.(curr_field).vbSPT.fastminAx50 = nanmedian(propdistr.(curr_field).vbSPT.fastminAx50);
    propdistr.(curr_field).vbSPT.fastminAx95 = data.vbSPT.fastminAx95(clusters_to_use);
    propPointEst.(curr_field).vbSPT.fastminAx95 = nanmedian(propdistr.(curr_field).vbSPT.fastminAx95);

    propdistr.(curr_field).vbSPT.slowbound50 = data.vbSPT.slowbound50(clusters_to_use);
    propdistr.(curr_field).vbSPT.slowbound95 = data.vbSPT.slowbound95(clusters_to_use);
    propdistr.(curr_field).vbSPT.fastbound50 = data.vbSPT.fastbound50(clusters_to_use);
    propdistr.(curr_field).vbSPT.fastbound95 = data.vbSPT.fastbound95(clusters_to_use);

    propdistr.(curr_field).vbSPT.slowCOM50 = data.vbSPT.slowCOM50(clusters_to_use,:);
    propdistr.(curr_field).vbSPT.slowCOM95 = data.vbSPT.slowCOM95(clusters_to_use,:);
    propdistr.(curr_field).vbSPT.fastCOM50 = data.vbSPT.fastCOM50(clusters_to_use,:);
    propdistr.(curr_field).vbSPT.fastCOM95 = data.vbSPT.fastCOM95(clusters_to_use,:);

    % lifetime
    propdistr.(curr_field).lifetime = data.lifetime(clusters_to_use);

    % recruitment intervals
    RI_subset = cell2mat(data.rec_int(clusters_to_use));
    propdistr.(curr_field).rec_int = RI_subset;

    % adhesion information
    propdistr.(curr_field).d2ad = data.adhesion.d2ad(clusters_to_use);
end
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
propPointEst.allLoc.KDE.area95 = nanmean(propdistr.allLoc.KDE.area95);

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

% lifetime
propdistr.lifetime = data.lifetime(clusters_to_use);

% % recruitment intervals
RI_subset = cell2mat(data.rec_int(clusters_to_use));
propdistr.rec_int = RI_subset;
end
