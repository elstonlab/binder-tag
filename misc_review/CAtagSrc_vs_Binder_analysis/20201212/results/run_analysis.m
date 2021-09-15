function run_analysis()
dataloc = 'inputdata/';
vbSPTDirLoc = 'inputdata/vbspt_result/';
vbSPTMetaDirLoc = 'inputdata/vbspt_source/';
dateprefix = '20201215';
parentFigDir = 'figures_20201215/';
nMinTracks = 10;
maxRIcutoff = 4;
numCPUs = 15;
analysis_mode = 'lvPALM';
cluster_segmentation(dataloc, dateprefix, parentFigDir, nMinTracks, ...
                     maxRIcutoff, vbSPTDirLoc, vbSPTMetaDirLoc,...
                     'temp_01_LLpardir', numCPUs, analysis_mode);
end

