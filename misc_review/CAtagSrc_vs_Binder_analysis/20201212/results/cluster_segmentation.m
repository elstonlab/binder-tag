% CLUSTER_SEGMENTATION.m
% Master controlling method for SPT cluster segmentation pipeline
% originally developed for the Tag-Src/Binder dataset. Intended to run
% on Longleaf for high-throughput.
%
% Workflow ==============================================================================================================
% 1. Obtain initial proposed clusters by Voronoi segmentation (Andronov L et al, Sci. Rep. 2016, 6:24084)
%    a. Simulate and Voronoi tessellate random distributions 50x per cell using cell ROI and known track abundance
%    b. Voronoi tessellate the experimental data using cell ROI and known track centroids
%    c. Use method described in paper to find polygon area threshold, fuse remaining polygons, and propose clusters
% 2. Spatial refinement (I): For all tracks un-assigned to clusters, assign them to a cluster if it diffuses into
%    the fused Voronoi polygon during its lifetime.
% 3. Temporal refinement (II): Determine recruitment intervals for all clusters; if there are any gaps > 4 seconds
%    (user-defined; 4s was empirically identified) then segregate the cluster into multiple pieces.
% **** The minimum # tracks in cluster can be set for steps 4 and 5.
%      Other filters were applied post-analysis to facilitate exploration of the filter effects. ****
% 4. Cluster characterization
% 5. Overall characterization
%
% 2018 May 08 / Mike Pablo ===============================================================================================

function cluster_segmentation(dataloc,dateprefix,parentFigDir,nMinTracks,...
                              maxRIcutoff,...
                              vbSPTDirLoc,vbSPTMetaDirLoc,...
                              longleafParDir,numcpus,analysis_mode)
% Input:
%   dataloc: string, location of lvPALM.mat and .roi files
%   dateprefix: string, recommended style '20180508_' for May 8 2018, etc.
%   parentFigDir: string, location to house all figures generated during analysis.
%                         recommended style 'figures_20180508/'
%   nMinTracks: double, number of minimum tracks required to register as a cluster during the
%                       initial Voronoi segmentation step.
%   maxRIcutoff: double, maximum time (in seconds) allowed between observing new molecules at a cluster
%                           before segregating a cluster
%   vbSPTDirLoc: string, location of vbSPT results
%   vbSPTMetaDirLoc: string, location of vbSPT metadata needed to link to raw data
%   longleafParDir: string, directory where parallel job information should be stored while running
%                   (critical for running multiple parallel jobs in the same folder on Longleaf)

assert(isstring(dataloc)||ischar(dataloc),'dataloc must be a string');
assert(isstring(dateprefix)||ischar(dateprefix),'dateprefix must be a string')
assert(isstring(parentFigDir)||ischar(parentFigDir),'parentFigDir must be a string')
assert(isnumeric(nMinTracks),'nMinTracks must be a numeric (default numeric type)')
assert(isnumeric(maxRIcutoff),'minRecruitmentGap must be a numeric (default numeric type)')
assert(isstring(vbSPTDirLoc)||ischar(vbSPTDirLoc),'vbSPTDirLoc must be a string')
assert(isstring(vbSPTMetaDirLoc)||ischar(vbSPTMetaDirLoc),'vbSPTMetaDirLoc must be a string')

dataloc = {dataloc}; % Encapsulate in cell array

% Suffixes of filenames to be generated/used
suffix.sim_den = 'sim_den_distr_50x';
suffix.exp_den = 'exp_den_distr';
suffix.thresh_file = 'voronoi_thresh';
suffix.voronoi_clusters = 'voronoi_clus';
suffix.refinementI_clusters = 'refinementI_clusTrackIDs';
suffix.refinementII_clusters = 'refinementII_clusTrackIDs';
suffix.cluster_props = 'measured_cluster_props';
suffix.cluster_analysis = strrep(sprintf('analyzed_cluster_props_maxRI%.2fs_minTrk%02d',maxRIcutoff,nMinTracks),'.','-');
suffix.overall_analysis = strrep(sprintf('analyzed_overall_props_maxRI%.2fs_minTrk%02d',maxRIcutoff,nMinTracks),'.','-');



% Set up dependencies -- there should be a check here for the lvPALM flag vs. the CDS flag.
if strcmpi(analysis_mode,'cds') % Adhesion datasets
    addpath('CDS_subroutines')
elseif strcmpi(analysis_mode,'lvpalm') % Non-adhesion datasets
    addpath('lvPALM_subroutines')
else
    error('invalid analysis mode set')
end
addpath(genpath('dependencies'));

% Set up parallel processing folder and settings
mkdir(longleafParDir);
pc=parcluster('local');
pc.JobStorageLocation = longleafParDir;
parpool(pc,numcpus);

mkdir(parentFigDir)

% Pipeline start
%simulate_uniform_density_distributions_50x(dataloc,dateprefix,suffix.sim_den); % outputs [dateprefix * suffix.sim_den] where * is tag or bin
%generate_density_distributions(dataloc,dateprefix,suffix.exp_den);             % outputs [dateprefix * suffix.exp_den] where * is tag or bin
%compare_voronoicell_size_distributions_sim50x(dateprefix,suffix.sim_den,suffix.exp_den,suffix.thresh_file,parentFigDir); % outputs [dateprefix suffix.thresh_file] as well as many figures
%generate_voronoiSegmented_clusters(dataloc,dateprefix,suffix.voronoi_clusters,suffix.thresh_file); % outputs [dateprefix * suffix.voronoi_clusters] where * is tag or bin
%cluster_refinement_I(dataloc,dateprefix,suffix.voronoi_clusters,suffix.refinementI_clusters); % outputs [dateprefix * suffix.refinementI_clusters] where * is tag or bin
%cluster_refinement_II(dataloc,dateprefix,suffix.refinementI_clusters,suffix.refinementII_clusters,maxRIcutoff,parentFigDir); % outputs [dateprefix suffix.refinementII_clusters] and some figures
%measure_cluster_props(dataloc,vbSPTDirLoc,vbSPTMetaDirLoc,dateprefix,suffix.refinementII_clusters,suffix.cluster_props); % outputs [dateprefix suffix.cluster_props]
analyze_cluster_props(dateprefix,suffix.cluster_props,suffix.cluster_analysis,nMinTracks,maxRIcutoff,parentFigDir); % outputs [dateprefix suffix.cluster_analysis] and lots of figures
analyze_overall_props(dataloc,vbSPTDirLoc,vbSPTMetaDirLoc,dateprefix,suffix.cluster_props,suffix.cluster_analysis,nMinTracks,suffix.overall_analysis,parentFigDir); % outputs [dateprefix suffix.overall_analysis] and lots of figures


end
