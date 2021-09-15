% CLUSTER_REFINEMENT_II
% Loads the refined cluster track IDs produced by cluster_refinement_I and
% the source datafiles to calculate recruitment intervals (time between
% addition of each new track to the cluster). We produce a plot of the current
% recruitment intervals and report the % of intervals below that value.
% Empirically, we refine the assignment based on the observation that
% ~90-95% of cluster recruitment % intervals are < 4 seconds, and the
% distribution strongly tailed off.
%
% Input:
%   datafolders: cell array of strings for all locations of raw data (*lvPALM.mat)
%   dateprefix and suffix_refinementII_clusters: strings for file naming, e.g.
%   '20180508_'    'clus_refinement_II'
%
%   suffix_refinementI_clusters: strings for file input, e.g.
%   'clus_refinement_I'
%
%   RIcutoff: numeric value for temporal cluster segmentation [seconds]
%
%   parentFigDir: where to save figures from this script
%
% Saves a .mat file to the working directory:
%   [dateprefix suffix_refinementII_clusters '.mat']
%
% Part of the cluster_segmentation.m pipeline.
% 2018 May 9 / Mike Pablo

function cluster_refinement_II(dataloc,dateprefix,suffix_refinementI_clusters,suffix_refinementII_clusters,RIcutoff,parentFigDir)
[all_path,all_file] = get_file_structures(dataloc);
ncells = numel(all_path.tag);

tagclusID = load([dateprefix 'tag' suffix_refinementI_clusters]);
binclusID = load([dateprefix 'bin' suffix_refinementI_clusters]);

figSubDir = 'RIs_pre_refinementII/';
parentFigDir = [parentFigDir figSubDir];
warning('off');
mkdir(parentFigDir);
warning('on');

tag_refined_cluster_track_IDs = cell(ncells,1);
bin_refined_cluster_track_IDs = cell(ncells,1);
parfor i=1:ncells
    curr_expdatafile_tag = [all_path.tag{i} '/' all_file.tag{i}];
    curr_expdatafile_bin = [all_path.bin{i} '/' all_file.bin{i}];

    tag_refined_cluster_track_IDs{i} = refine_cluster_results(tagclusID.cluster_track_IDs{i},...
                                                              curr_expdatafile_tag,RIcutoff,...
                                                              parentFigDir,sprintf('refinementII_RI_tag_%02d',i));
    bin_refined_cluster_track_IDs{i} = refine_cluster_results(binclusID.cluster_track_IDs{i},...
                                                              curr_expdatafile_bin,RIcutoff,...
                                                              parentFigDir,sprintf('refinementII_RI_bin_%02d',i));
end
save([dateprefix suffix_refinementII_clusters],'tag_refined_cluster_track_IDs','bin_refined_cluster_track_IDs');

end

function cluster_track_IDs = refine_cluster_results(clusTrackIDs,expdatafile,RIcutoff,parentFigDir,figName)
% Calculates recruitment interval data and reassigns the clusters
% If recruitment interval is greater than RIcutoff, then consider it a new cluster

fprintf('Refining clustering for:\n\t%s\n',expdatafile)
expdata = load(expdatafile,'smLinked');
initial_nclus = numel(clusTrackIDs);
cluster_track_IDs = {};

original_RI_distribution = cell(numel(clusTrackIDs),1);

%% Reassign tracks to clusters.
for j=1:initial_nclus
    trackIDs = clusTrackIDs{j};
    track_starts = cellfun(@(x) x(1,end),expdata.smLinked(trackIDs));

    % Determine track starts and corresponding IDs for the cluster.
    trackStartandID = [track_starts(:), trackIDs(:)]; % trackIDs is 1xN vector, starts is Nx1
    trackStartandID=sortrows(trackStartandID,1); % Sort by start point
    all_recruitment_intervals = (trackStartandID(2:end,1)-trackStartandID(1:end-1,1))*0.02;
    original_RI_distribution{j} = all_recruitment_intervals;
    % Does this cluster need temporal separation?
    if any(all_recruitment_intervals>RIcutoff)
        temp_RIs = [all_recruitment_intervals;nan];
        subcluster_assignments = nan(numel(trackIDs),1);
        curr_subcluster = 1;

        % Assign all the tracks into subclusters.
        for k=1:(numel(temp_RIs)-1) % last entry in tempRIs is edge case
           if temp_RIs(k) < RIcutoff
               subcluster_assignments(k) = curr_subcluster;
           else
               subcluster_assignments(k) = curr_subcluster;
               curr_subcluster = curr_subcluster+1;
           end
        end
        subcluster_assignments(end) = curr_subcluster;
        nsubclusters = curr_subcluster;

        for k=1:nsubclusters
            subcluster_idxs  = subcluster_assignments==k;
            cluster_track_IDs = [cluster_track_IDs,trackStartandID(subcluster_idxs,2)];
        end
    else % the whole cluster was 'normal', we don't even need to use re-sorted information.
        cluster_track_IDs = [cluster_track_IDs,trackIDs];
    end
end

original_RI_distribution = cell2mat(original_RI_distribution);

% Plot the original RI distribution
f=figure('position',[1000, 1007, 608, 331]);
histogram(original_RI_distribution,'binedges',0:0.02:max(original_RI_distribution))
xlabel('recruitment interval (s)')
ylabel('count')
percent_below_threshold = (sum(original_RI_distribution<RIcutoff)/numel(original_RI_distribution))*100;
title(sprintf('%.2f %% of raw RIs < %.2f',percent_below_threshold,RIcutoff));
set(gca,'fontsize',14,'fontname','arial','tickdir','out')

savename = [parentFigDir figName];
savefig(f,[savename '.fig']);
print(f,'-dtiff',[savename '.tif'],'-r300');

set(gca,'xlim',[0 RIcutoff]);
savename = [parentFigDir figName];
savefig(f,[savename '_zoom.fig']);
print(f,'-dtiff',[savename '_zoom.tif'],'-r300');

close(f);
fprintf('For %s, %.2f %% of input RIs < %.2f\n',figName,percent_below_threshold,RIcutoff);
end


function [all_path,all_file] = get_file_structures(datafolders)
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
end

function [vbSPTfile,vbSPTmetafile] = get_vbSPT_file_structures(ncells)
vbSPTfile.tag = cell(ncells,1);
vbSPTfile.bin = cell(ncells,1);
vbSPTmetafile.tag = cell(ncells,1);
vbSPTmetafile.bin = cell(ncells,1);

vbSPTdata_dir = './vbspt/results/';
vbSPTmeta_dir = './vbspt/vbspt_source/';
for i=1:ncells
    vbSPTfile.tag{i} = sprintf('%stag_%02d_HMManalysis_hidden2.mat',vbSPTdata_dir,i);
    vbSPTfile.bin{i} = sprintf('%sbin_%02d_HMManalysis_hidden2.mat',vbSPTdata_dir,i);
    vbSPTmetafile.tag{i} = sprintf('%stag_%02d.mat',vbSPTmeta_dir,i);
    vbSPTmetafile.bin{i} = sprintf('%sbin_%02d.mat',vbSPTmeta_dir,i);
end


end
