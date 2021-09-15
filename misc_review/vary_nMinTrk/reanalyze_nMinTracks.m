% See explanation.txt

function reanalyze_nMinTracks()

nMinTracksToTest = [3:20];

dependency_paths = {'C:\Users\mike\Documents\GitHub\binder-tag\initial_review\CAtagSrc_vs_Binder_analysis\20201212\results\lvPALM_subroutines\',...
                    'C:\Users\mike\Documents\GitHub\binder-tag\initial_review\CAtagSrc_vs_Binder_analysis\20201212\results\dependencies\',...
                    'J:\Elston Lab\MikePablo\20180504-bindertag\20190410_rerun_pipeline_and_binomanalysis\result'};
initialize_dependencies(dependency_paths);


basesources = {'J:\Elston Lab\MikePablo\20180716 - satadrive backup\20180524_bintag_0328-0405_mergecopy\',...
               'J:\Elston Lab\MikePablo\20180716 - satadrive backup\20180628_more_noAd\to_LL\01092018_100pg\',...
               'J:\Elston Lab\MikePablo\20180716 - satadrive backup\20180628_more_noAd\to_LL\01282018_100pg\',...
               'J:\Elston Lab\MikePablo\20180716 - satadrive backup\20180628_more_noAd\to_LL\06172018_100pg\'};
dataloc = cell(4,1);
vbSPTDirLoc = cell(4,1);
vbSPTMetaDirLoc = cell(4,1);
for i=1:4
   dataloc{i} = [basesources{i}, 'inputdata\'];
   vbSPTDirLoc{i} = [basesources{i}, 'inputdata\vbspt\results\'];
   vbSPTMetaDirLoc{i} = [basesources{i}, 'inputdata\vbspt\vbspt_source\'];
end
           
% dataloc = 'C:\Users\mike\Documents\GitHub\binder-tag\initial_review\CAtagSrc_vs_Binder_analysis\20201212\results\inputdata\';
% vbSPTDirLoc = 'C:\Users\mike\Documents\GitHub\binder-tag\initial_review\CAtagSrc_vs_Binder_analysis\20201212\results\inputdata\vbspt_result\';
% vbSPTMetaDirLoc = 'C:\Users\mike\Documents\GitHub\binder-tag\initial_review\CAtagSrc_vs_Binder_analysis\20201212\results\inputdata\vbspt_source\';
longleafParDir = 'none';


%0328-0405 - inputdata\vbspt\vbspt_source
%       inputdata\vbspt\results
%0109 - inputdata\vbspt\vbspt_source
%       inputdata\vbspt\results
%0128 - inputdata\vbspt\vbspt_source
%       inputdata\vbspt\results
%0617 - inputdata\vbspt\vbspt_source
%       inputdata\vbspt\results

maxRIcutoff = 4;
numcpus = 1;
analysis_mode = 'lvPALM';

dateprefix = {'20180524_0328-0405_noAd_100pg_',...
              '20180628_noAd_0109_100pg',...
              '20180628_noAd_0128_100pg',...
              '20180628_noAd_0617_100pg',...
              };
for j=1:4
    for i=1:numel(nMinTracksToTest)
       nMinTracks = nMinTracksToTest(i);
       
       
    if ~(j==1 && nMinTracks >= 3 && nMinTracks <=6)
        continue;
    end
       nMinTrackDir = sprintf('nMinTrack_%02d', nMinTracks);
       parentFigDir = 'none/';
       finish_cluster_segmentation(dataloc{j},dateprefix{j},parentFigDir,nMinTracks,...
                                   maxRIcutoff,...
                                   vbSPTDirLoc{j},vbSPTMetaDirLoc{j},...
                                   longleafParDir,numcpus,analysis_mode,...
                                   nMinTrackDir);
    end
end

cleanup_dependencies(dependency_paths);
end

function initialize_dependencies(dependency_paths)

n_paths = numel(dependency_paths);
for i=1:n_paths
   addpath(genpath(dependency_paths{i})); 
end
end

function cleanup_dependencies(dependency_paths)
n_paths = numel(dependency_paths);
for i=1:n_paths
   rmpath(genpath(dependency_paths{i})); 
end
end

function finish_cluster_segmentation(dataloc,dateprefix,parentFigDir,nMinTracks,...
                                     maxRIcutoff,...
                                     vbSPTDirLoc,vbSPTMetaDirLoc,...
                                     longleafParDir,numcpus,analysis_mode,...
                                     nMinTrackDir)

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
suffix.cluster_analysis = strrep(sprintf('analyzed_cluster_props_maxRI%.2fs_minTrk%02d_maxClusSize-0-60um_SPTacqTime60',maxRIcutoff,nMinTracks),'.','-');
suffix.overall_analysis = strrep(sprintf('analyzed_overall_props_maxRI%.2fs_minTrk%02d_maxClusSize-0-60um_SPTacqTime60',maxRIcutoff,nMinTracks),'.','-');

mkdir(nMinTrackDir);
cd(nMinTrackDir);

fprintf('Working in %s on %s\n', nMinTrackDir, dataloc{1});

mkdir(parentFigDir);
analyze_cluster_props(dateprefix,suffix.cluster_props,suffix.cluster_analysis,nMinTracks,maxRIcutoff,parentFigDir); % outputs [dateprefix suffix.cluster_analysis] and lots of figures
analyze_overall_props(dataloc,vbSPTDirLoc,vbSPTMetaDirLoc,dateprefix,suffix.cluster_props,suffix.cluster_analysis,nMinTracks,suffix.overall_analysis,parentFigDir); % outputs [dateprefix suffix.overall_analysis] and lots of figures
 
cd ..
end