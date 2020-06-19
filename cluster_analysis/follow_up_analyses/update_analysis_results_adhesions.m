function update_analysis_results_adhesions()
% Organize all the output results from cluster analysis at this directory.
% This code requires tha the cluster analysis code was already performed.
sourcedir = '/Volumes/SATAdrive/20190521_pxn_rerun/';

% We had one sets of separately-conducted analyses to correct.
% There are the "measured cluster props" results and "analyzed cluster props" results:
meas1=load([sourcedir '20190520_measured_cluster_props.mat']);
allmeas.tag_prop=[meas1.tag_prop];
allmeas.bin_prop=[meas1.bin_prop];

anal1=load([sourcedir '20190520_analyzed_cluster_props_maxRI4-00s_minTrk10_maxClusSize0-60um_SPTacqTime10.mat']);
allanal_propdists_tag=[anal1.propdists_byAd.tag];
allanal_propdists_bin=[anal1.propdists_byAd.bin];

% We needed to re-load the experimental data for these follow-up analyses
[path1_tag,file1_tag]=batchGetPath('/Volumes/med/Elston Lab/MikePablo/20180716 - satadrive backup/binder_tag_datacopy/spt_bt_datacopy/pxn_100pgdox_cds_removedBadCells/','RightCali_lvPALM_corrDataStruct','mat');
[path1_bin,file1_bin]=batchGetPath('/Volumes/med/Elston Lab/MikePablo/20180716 - satadrive backup/binder_tag_datacopy/spt_bt_datacopy/pxn_100pgdox_cds_removedBadCells/','Left_lvPALM_corrDataStruct','mat');

allpath_tag = [path1_tag];
allpath_bin = [path1_bin];

allfile_tag = [file1_tag];
allfile_bin = [file1_bin];

ncells = numel(allpath_tag);
% The analysis was conducted in parallel, but the parfor could be converted into a simple for loop.
parfor i=1:ncells
    fprintf('cell %i\n',i)
    curr_tagexpdatafile = fullfile(allpath_tag{i},allfile_tag{i});
    curr_binexpdatafile = fullfile(allpath_bin{i},allfile_bin{i});

    [allanal_propdists_tag{i},allanal_propdists_bin{i}] = annotate_with_filter_info(allanal_propdists_tag{i},allanal_propdists_bin{i},...
                                                                                    allmeas.tag_prop{i},allmeas.bin_prop{i},...
                                                                                    curr_tagexpdatafile,curr_binexpdatafile);



end
% We perform this renaming just because indexing within the fields of a struct is problematic in a parfor loop
allanal.propdists_byAd.tag = allanal_propdists_tag;
allanal.propdists_byAd.bin = allanal_propdists_bin;
save('annotated_propdists','allanal');
end

function [tagpropdistInst,binpropdistInst] = annotate_with_filter_info(tagpropdistInst,binpropdistInst,tagmeasInst,binmeasInst,tagexpdatafile,binexpdatafile)
% We're going to annotate the data with additional information, relevant for filtering
% The main ones are as follows:
% clus_lifetime_vector    :     non-zero if a track was present in the cluster
%                               during the frames. Has length equal to the # frames in the cluster lifetime.
%                               Can be >1 if multiple tracks were observed at the same time.
% CLVgt2                  : Stands for 'clus_lifetime_vector greater than 2'; checks if two or more tracks
%                           were observed at the same. Boolean.
% startAfter4Sec          : Checks if the cluster started 4 seconds after imaging started. Boolean.


expdata_tag = load(tagexpdatafile);
expdata_bin = load(binexpdatafile);

adfields = {'NA','FC','FA','off','on','small','large'};

tagclusIDs.NA = (tagmeasInst.adhesion.d2ad < 0.1067*2) & tagmeasInst.N >= 10 & (tagmeasInst.adsize'>= 5) & (tagmeasInst.adsize'<22);
tagclusIDs.FC = (tagmeasInst.adhesion.d2ad < 0.1067*2) & tagmeasInst.N >= 10 & (tagmeasInst.adsize'>=22) & (tagmeasInst.adsize'<88);
tagclusIDs.FA = (tagmeasInst.adhesion.d2ad < 0.1067*2) & tagmeasInst.N >= 10 & (tagmeasInst.adsize'>=88) & (tagmeasInst.adsize'<440);
tagclusIDs.off = ~(tagmeasInst.adhesion.d2ad < 0.1067*2) & tagmeasInst.N >= 10;
tagclusIDs.on = tagclusIDs.NA | tagclusIDs.FC | tagclusIDs.FA;
tagclusIDs.small = (tagmeasInst.adhesion.d2ad < 0.1067*2) & tagmeasInst.N >= 10 & (tagmeasInst.adsize'>= 5) & (tagmeasInst.adsize'<88);
tagclusIDs.large = (tagmeasInst.adhesion.d2ad < 0.1067*2) & tagmeasInst.N >= 10 & (tagmeasInst.adsize'>= 88) & (tagmeasInst.adsize'<440);
binclusIDs.NA = (binmeasInst.adhesion.d2ad < 0.1067*2) & binmeasInst.N >= 10 & (binmeasInst.adsize'>= 5) & (binmeasInst.adsize'<22);
binclusIDs.FC = (binmeasInst.adhesion.d2ad < 0.1067*2) & binmeasInst.N >= 10 & (binmeasInst.adsize'>=22) & (binmeasInst.adsize'<88);
binclusIDs.FA = (binmeasInst.adhesion.d2ad < 0.1067*2) & binmeasInst.N >= 10 & (binmeasInst.adsize'>=88) & (binmeasInst.adsize'<440);
binclusIDs.off = ~(binmeasInst.adhesion.d2ad < 0.1067*2) & binmeasInst.N >= 10;
binclusIDs.on = binclusIDs.NA | binclusIDs.FC | binclusIDs.FA;
binclusIDs.small = (binmeasInst.adhesion.d2ad < 0.1067*2) & binmeasInst.N >= 10 & (binmeasInst.adsize'>= 5) & (binmeasInst.adsize'<88);
binclusIDs.large = (binmeasInst.adhesion.d2ad < 0.1067*2) & binmeasInst.N >= 10 & (binmeasInst.adsize'>= 88) & (binmeasInst.adsize'<440);


for adfieldid = 1:numel(adfields)
    curradfield = adfields{adfieldid};



    tagmeasInst.clus_to_use = tagclusIDs.(curradfield);
    binmeasInst.clus_to_use = binclusIDs.(curradfield);

    tagClusTrackIDs = tagmeasInst.trackIDs(tagmeasInst.clus_to_use);
    binClusTrackIDs = binmeasInst.trackIDs(binmeasInst.clus_to_use);

    tagClusAdFrames = tagmeasInst.adhesion.adframe(tagmeasInst.clus_to_use);
    binClusAdFrames = binmeasInst.adhesion.adframe(binmeasInst.clus_to_use);

    assert(numel(tagClusTrackIDs) == numel(tagpropdistInst.(curradfield).N),'mismatch # tag clusters');
    assert(numel(binClusTrackIDs) == numel(binpropdistInst.(curradfield).N),'mismatch # bin clusters');

    nclus_tag = numel(tagClusTrackIDs);
    nclus_bin = numel(binClusTrackIDs);

    tagpropdistInst.(curradfield).clus_lifetime_vector = cell(nclus_tag,1);
    binpropdistInst.(curradfield).clus_lifetime_vector = cell(nclus_bin,1);
    tagpropdistInst.(curradfield).clusAdFrames = tagClusAdFrames;
    binpropdistInst.(curradfield).clusAdFrames = binClusAdFrames;
    for i=1:nclus_tag
        % Get all the tracks
        curr_tag_track_ids = tagClusTrackIDs{i};

        nTag = numel(curr_tag_track_ids);
        curr_tag_firstframes = cellfun(@(x) x.allfeature(1,end),expdata_tag.corrDataStruct(curr_tag_track_ids));
        curr_tag_lastframes = cellfun(@(x) x.allfeature(end,end),expdata_tag.corrDataStruct(curr_tag_track_ids));

        clus_first_frame =  min(curr_tag_firstframes);
        clus_last_frame =  max(curr_tag_lastframes);

        clus_lifetime_vector = zeros(numel(clus_first_frame:clus_last_frame),1);
        for j=1:nTag
            adjusted_tag_track_start = (curr_tag_firstframes(j)-clus_first_frame+1);
            adjusted_tag_track_end = (curr_tag_lastframes(j)-clus_first_frame+1);

            clus_lifetime_vector( adjusted_tag_track_start : adjusted_tag_track_end) = clus_lifetime_vector( adjusted_tag_track_start : adjusted_tag_track_end) + 1;
        end
        tagpropdistInst.(curradfield).clus_lifetime_vector{i} = clus_lifetime_vector;
        tagpropdistInst.(curradfield).CLVgt2(i) = any(clus_lifetime_vector>2);
        tagpropdistInst.(curradfield).startAfter4Sec = all(adjusted_tag_track_start>(4/.02));
    end
    for i=1:nclus_bin
        % Get all the tracks
        curr_bin_track_ids = binClusTrackIDs{i};

        nBin = numel(curr_bin_track_ids);
        curr_bin_firstframes = cellfun(@(x) x.allfeature(1,end),expdata_bin.corrDataStruct(curr_bin_track_ids));
        curr_bin_lastframes = cellfun(@(x) x.allfeature(end,end),expdata_bin.corrDataStruct(curr_bin_track_ids));

        clus_first_frame =  min(curr_bin_firstframes);
        clus_last_frame =  max(curr_bin_lastframes);

        clus_lifetime_vector = zeros(numel(clus_first_frame:clus_last_frame),1);
        for j=1:nBin
            adjusted_bin_track_start = (curr_bin_firstframes(j)-clus_first_frame+1);
            adjusted_bin_track_end = (curr_bin_lastframes(j)-clus_first_frame+1);

            clus_lifetime_vector( adjusted_bin_track_start : adjusted_bin_track_end) = clus_lifetime_vector( adjusted_bin_track_start : adjusted_bin_track_end) + 1;
        end

        binpropdistInst.(curradfield).clus_lifetime_vector{i} = clus_lifetime_vector;
        binpropdistInst.(curradfield).CLVgt2(i) = any(clus_lifetime_vector>2);
        binpropdistInst.(curradfield).startAfter4Sec = all(adjusted_bin_track_start>(4/.02));
    end

end
end
