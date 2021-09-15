function redo_postprocessing_LL()
nMinTracksToTest = [10];
sourcedir = 'data/';
addpath(genpath('dependencies'));
addpath(genpath('lvPALM_subroutines'));
inputdatapaths = {'data/20180524_bintag_0328-0405_mergecopy/inputdata/',...
                  'data/01092018_100pg/inputdata/',...
                  'data/01282018_100pg/inputdata/',...
                  'data/06172018_100pg/inputdata/'};

               
for i=1:numel(nMinTracksToTest)
    nMinTracks = nMinTracksToTest(i);
    nMinTrackDir = sprintf('nMinTrack_%02d/', nMinTracks);
    resultDir = sprintf('nMinTrack_result_%02d/', nMinTracks);
    analfile1 = [nMinTrackDir, sprintf('20180524_0328-0405_noAd_100pg_analyzed_cluster_props_maxRI4-00s_minTrk%02d_maxClusSize-0-60um_SPTacqTime60.mat', nMinTracks)];
    analfile2 = [nMinTrackDir, sprintf('20180628_noAd_0109_100pganalyzed_cluster_props_maxRI4-00s_minTrk%02d_maxClusSize-0-60um_SPTacqTime60.mat', nMinTracks)];
    analfile3 = [nMinTrackDir, sprintf('20180628_noAd_0128_100pganalyzed_cluster_props_maxRI4-00s_minTrk%02d_maxClusSize-0-60um_SPTacqTime60.mat', nMinTracks)];
    analfile4 = [nMinTrackDir, sprintf('20180628_noAd_0617_100pganalyzed_cluster_props_maxRI4-00s_minTrk%02d_maxClusSize-0-60um_SPTacqTime60.mat', nMinTracks)];

    measfile1 = [sourcedir '20180524_0328-0405_noAd_100pg_measured_cluster_props.mat'];
    measfile2 = [sourcedir '20180628_noAd_0109_100pgmeasured_cluster_props.mat'];
    measfile3 = [sourcedir '20180628_noAd_0128_100pgmeasured_cluster_props.mat'];
    measfile4 = [sourcedir '20180628_noAd_0617_100pgmeasured_cluster_props.mat'];

    analfiles = {analfile1, analfile2, analfile3, analfile4};
    measfiles = {measfile1, measfile2, measfile3, measfile4};
    
    savename = [resultDir, sprintf('annotated_propdists_wRecToClus_minTrk%02d',nMinTracks)];
    
    update_analysis_results(measfiles,...
                            analfiles,...
                            inputdatapaths,...
                            savename,...
                            nMinTracks)
end
               
end

function update_analysis_results(measfiles,...
                                 analfiles,...
                                 inputdatapaths,...
                                 savename,...
                                 nMinTracks)
                           
% Organize all the output results from cluster analysis at this directory.
% This code requires tha the cluster analysis code was already performed.


% We had four sets of separately-conducted analyses to correct.
% There are the "measured cluster props" results and "analyzed cluster props" results:
meas1=load(measfiles{1});
meas2=load(measfiles{2});
meas3=load(measfiles{3});
meas4=load(measfiles{4});
allmeas.tag_prop=[meas1.tag_prop;meas2.tag_prop;meas3.tag_prop;meas4.tag_prop];
allmeas.bin_prop=[meas1.bin_prop;meas2.bin_prop;meas3.bin_prop;meas4.bin_prop];

anal1=load(analfiles{1});
anal2=load(analfiles{2});
anal3=load(analfiles{3});
anal4=load(analfiles{4});
allanal_propdists_tag=[anal1.propdists.tag;anal2.propdists.tag;anal3.propdists.tag;anal4.propdists.tag];
allanal_propdists_bin=[anal1.propdists.bin;anal2.propdists.bin;anal3.propdists.bin;anal4.propdists.bin];

% We needed to re-load the experimental data for these follow-up analyses
[path1_tag,file1_tag]=batchGetPath(inputdatapaths{1},'RightCali_lvPALM','mat');
[path2_tag,file2_tag]=batchGetPath(inputdatapaths{2},'RightCali_lvPALM','mat');
[path3_tag,file3_tag]=batchGetPath(inputdatapaths{3},'RightCali_lvPALM','mat');
[path4_tag,file4_tag]=batchGetPath(inputdatapaths{4},'RightCali_lvPALM','mat');

[path1_bin,file1_bin]=batchGetPath(inputdatapaths{1},'Left_lvPALM','mat');
[path2_bin,file2_bin]=batchGetPath(inputdatapaths{2},'Left_lvPALM','mat');
[path3_bin,file3_bin]=batchGetPath(inputdatapaths{3},'Left_lvPALM','mat');
[path4_bin,file4_bin]=batchGetPath(inputdatapaths{4},'Left_lvPALM','mat');

% On Windows, '\' is used for paths, but on Unix/Mac, '/' is used.
% When sorting string arrays, they do not have same priority
% (e.g. '/' is higher priority than numeric characters, but
%       '\' is lower priority than numeric characters)
% This can result in mismatched files; we can correct it by sorting on reformatted strings
[~,idx]=sort(cellfun(@(x) strrep(x,'\','/'), path1_tag,'uniformoutput',false));
path1_tag = path1_tag(idx);
file1_tag = file1_tag(idx);
[~,idx]=sort(cellfun(@(x) strrep(x,'\','/'), path2_tag,'uniformoutput',false));
path2_tag = path2_tag(idx);
file2_tag = file2_tag(idx);
[~,idx]=sort(cellfun(@(x) strrep(x,'\','/'), path3_tag,'uniformoutput',false));
path3_tag = path3_tag(idx);
file3_tag = file3_tag(idx);
[~,idx]=sort(cellfun(@(x) strrep(x,'\','/'), path4_tag,'uniformoutput',false));
path4_tag = path4_tag(idx);
file4_tag = file4_tag(idx);
[~,idx]=sort(cellfun(@(x) strrep(x,'\','/'), path1_bin,'uniformoutput',false));
path1_bin = path1_bin(idx);
file1_bin = file1_bin(idx);
[~,idx]=sort(cellfun(@(x) strrep(x,'\','/'), path2_bin,'uniformoutput',false));
path2_bin = path2_bin(idx);
file2_bin = file2_bin(idx);
[~,idx]=sort(cellfun(@(x) strrep(x,'\','/'), path3_bin,'uniformoutput',false));
path3_bin = path3_bin(idx);
file3_bin = file3_bin(idx);
[~,idx]=sort(cellfun(@(x) strrep(x,'\','/'), path4_bin,'uniformoutput',false));
path4_bin = path4_bin(idx);
file4_bin = file4_bin(idx);

allpath_tag = [path1_tag;path2_tag;path3_tag;path4_tag];
allpath_bin = [path1_bin;path2_bin;path3_bin;path4_bin];

allfile_tag = [file1_tag;file2_tag;file3_tag;file4_tag];
allfile_bin = [file1_bin;file2_bin;file3_bin;file4_bin];

ncells = numel(allpath_tag);
% The analysis was conducted in parallel, but the parfor could be converted into a simple for loop.
% This analysis seems extremely slow on my computer; i think I run out of memory?
% Maybe if i rewrite it to only load parts of the data at a time..
% alternatively i h ave to do this via longleaf.
for i=1:ncells
    fprintf('cell %i\n',i)
    curr_tagexpdatafile = fullfile(allpath_tag{i},allfile_tag{i});
    curr_binexpdatafile = fullfile(allpath_bin{i},allfile_bin{i});

    [allanal_propdists_tag{i},allanal_propdists_bin{i}] = annotate_with_filter_info(allanal_propdists_tag{i},allanal_propdists_bin{i},...
                                                                                    allmeas.tag_prop{i},allmeas.bin_prop{i},...
                                                                                    curr_tagexpdatafile,curr_binexpdatafile,...
                                                                                    nMinTracks);



end
% We perform this renaming just because indexing within the fields of a struct is problematic in a parfor loop
allanal.propdists.tag = allanal_propdists_tag;
allanal.propdists.bin = allanal_propdists_bin;

%savename
save(savename,'allanal');
end

function [tagpropdistInst,binpropdistInst] = annotate_with_filter_info(tagpropdistInst,binpropdistInst,tagmeasInst,binmeasInst,tagexpdatafile,binexpdatafile,nMinTracks)
% We're going to annotate the data with additional informatino, relevant for filtering
% This includes the following:
% clus_lifetime_vector    :     non-zero if a track was present in the cluster
%                               during the frames. Has length equal to the # frames in the cluster lifetime.
%                               Can be >1 if multiple tracks were observed at the same time.
% CLVgt2                  : Stands for 'clus_lifetime_vector greater than 2'; checks if two or more tracks
%                           were observed at the same. Boolean.
% startAfter4Sec          : Checks if the cluster started 4 seconds after imaging started. Boolean.
% trackCOM                : Track-basec center of mass. Another way to estimate the cluster's position.

% We also include additional measurements for the recToClus structure,
% reporting the absolute numbers of moleculse in the cluster and the absolute number starting in the central region.
% These are important for estimating a simulated random level of recruitment.

% For the paper, we required clusters to have CLVgt2 == false and startAfter4Sec == true
% in addition to the other filters (e.g. >=10 tracks, spatial/temporal refinement)

expdata_tag = load(tagexpdatafile);
expdata_bin = load(binexpdatafile);

%tagmeasInst.clus_to_use = tagmeasInst.N >= 10;
%binmeasInst.clus_to_use = binmeasInst.N >= 10;
tagmeasInst.clus_to_use = tagmeasInst.N >= nMinTracks;
binmeasInst.clus_to_use = binmeasInst.N >= nMinTracks;



tagClusTrackIDs = tagmeasInst.trackIDs(tagmeasInst.clus_to_use);
binClusTrackIDs = binmeasInst.trackIDs(binmeasInst.clus_to_use);
tagClusTrackBounds = tagmeasInst.allLoc.KDEbound50(tagmeasInst.clus_to_use);
binClusTrackBounds = binmeasInst.allLoc.KDEbound50(binmeasInst.clus_to_use);

assert(numel(tagClusTrackIDs) == numel(tagpropdistInst.N),'mismatch # tag clusters');
assert(numel(binClusTrackIDs) == numel(binpropdistInst.N),'mismatch # bin clusters');

nclus_tag = numel(tagClusTrackIDs);
nclus_bin = numel(binClusTrackIDs);

tagpropdistInst.clus_lifetime_vector = cell(nclus_tag,1);
binpropdistInst.clus_lifetime_vector = cell(nclus_bin,1);
for i=1:nclus_tag
    % Get all the tracks
    curr_tag_track_ids = tagClusTrackIDs{i};

    nTag = numel(curr_tag_track_ids);
    curr_tag_firstframes = cellfun(@(x) x(1,end),expdata_tag.smLinked(curr_tag_track_ids));
    curr_tag_lastframes = cellfun(@(x) x(end,end),expdata_tag.smLinked(curr_tag_track_ids));

    clus_first_frame =  min(curr_tag_firstframes);
    clus_last_frame =  max(curr_tag_lastframes);

    clus_lifetime_vector = zeros(numel(clus_first_frame:clus_last_frame),1);
    for j=1:nTag
        adjusted_tag_track_start = (curr_tag_firstframes(j)-clus_first_frame+1);
        adjusted_tag_track_end = (curr_tag_lastframes(j)-clus_first_frame+1);

        clus_lifetime_vector( adjusted_tag_track_start : adjusted_tag_track_end) = clus_lifetime_vector( adjusted_tag_track_start : adjusted_tag_track_end) + 1;
    end
    tagpropdistInst.clus_lifetime_vector{i} = clus_lifetime_vector;
    tagpropdistInst.CLVgt2(i) = any(clus_lifetime_vector>2);
    tagpropdistInst.startAfter4Sec(i) = all(adjusted_tag_track_start>(4/.02));
    tagpropdistInst.trackCOM(i,:) = [nanmean(cell2mat(cellfun(@(x) x(:,1),expdata_tag.smLinked(curr_tag_track_ids),'uniformoutput',false))) , ...
                                nanmean(cell2mat(cellfun(@(x) x(:,2),expdata_tag.smLinked(curr_tag_track_ids),'uniformoutput',false)))];

    % recrutiment to cluster
    if all(~isnan(tagClusTrackBounds{i}(:)))
        warning('off')
        cluster_region = polyshape(tagClusTrackBounds{i});
        warning('on')

        first_coords = cell2mat(cellfun(@(x) x(1,1:2),expdata_tag.smLinked(curr_tag_track_ids),'uniformoutput',false))*0.1067;

        is_rec_to_clus = isinterior(cluster_region,first_coords);
        KDEr50_nRecToClus_maskBound = sum(is_rec_to_clus);
        KDEr50_nInClus_maskBound = nTag;
        KDEr50_fracRecToClus_maskBound = KDEr50_nRecToClus_maskBound/KDEr50_nInClus_maskBound;

        tagpropdistInst.recToClus.frac(i) = KDEr50_fracRecToClus_maskBound;
        tagpropdistInst.recToClus.N_recTracks(i) = KDEr50_nRecToClus_maskBound;
        tagpropdistInst.recToClus.N_allTracks(i) = KDEr50_nInClus_maskBound;
    else
        tagpropdistInst.recToClus.frac(i) = nan;
        tagpropdistInst.recToClus.N_recTracks(i) = nan;
        tagpropdistInst.recToClus.N_allTracks(i) = nTag;
    end
end
for i=1:nclus_bin
    % Get all the tracks
    curr_bin_track_ids = binClusTrackIDs{i};

    nBin = numel(curr_bin_track_ids);
    curr_bin_firstframes = cellfun(@(x) x(1,end),expdata_bin.smLinked(curr_bin_track_ids));
    curr_bin_lastframes = cellfun(@(x) x(end,end),expdata_bin.smLinked(curr_bin_track_ids));

    clus_first_frame =  min(curr_bin_firstframes);
    clus_last_frame =  max(curr_bin_lastframes);

    clus_lifetime_vector = zeros(numel(clus_first_frame:clus_last_frame),1);
    for j=1:nBin
        adjusted_bin_track_start = (curr_bin_firstframes(j)-clus_first_frame+1);
        adjusted_bin_track_end = (curr_bin_lastframes(j)-clus_first_frame+1);

        clus_lifetime_vector( adjusted_bin_track_start : adjusted_bin_track_end) = clus_lifetime_vector( adjusted_bin_track_start : adjusted_bin_track_end) + 1;
    end

    binpropdistInst.clus_lifetime_vector{i} = clus_lifetime_vector;
    binpropdistInst.CLVgt2(i) = any(clus_lifetime_vector>2);
    binpropdistInst.startAfter4Sec(i) = all(adjusted_bin_track_start>(4/.02));
    binpropdistInst.trackCOM(i,:) = [nanmean(cell2mat(cellfun(@(x) x(:,1),expdata_bin.smLinked(curr_bin_track_ids),'uniformoutput',false))) , ...
                                nanmean(cell2mat(cellfun(@(x) x(:,2),expdata_bin.smLinked(curr_bin_track_ids),'uniformoutput',false)))];


    % recrutiment to cluster
    if all(~isnan(binClusTrackBounds{i}(:)))
        warning('off')
        cluster_region = polyshape(binClusTrackBounds{i});
        warning('on')

        first_coords = cell2mat(cellfun(@(x) x(1,1:2),expdata_bin.smLinked(curr_bin_track_ids),'uniformoutput',false))*0.1067;

        is_rec_to_clus = isinterior(cluster_region,first_coords);
        KDEr50_nRecToClus_maskBound = sum(is_rec_to_clus);
        KDEr50_nInClus_maskBound = nBin;
        KDEr50_fracRecToClus_maskBound = KDEr50_nRecToClus_maskBound/KDEr50_nInClus_maskBound;

        binpropdistInst.recToClus.frac(i) = KDEr50_fracRecToClus_maskBound;
        binpropdistInst.recToClus.N_recTracks(i) = KDEr50_nRecToClus_maskBound;
        binpropdistInst.recToClus.N_allTracks(i) = KDEr50_nInClus_maskBound;
    else
        binpropdistInst.recToClus.frac(i) = nan;
        binpropdistInst.recToClus.N_recTracks(i) = nan;
        binpropdistInst.recToClus.N_allTracks(i) = nBin;
    end
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
%

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
