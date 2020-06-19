# cluster_analysis

## Overview
The overall workflow is as follows. Steps 1-3 are for setting up inputs to the pipeline.
1. Single particle tracking of tagSrc and Binder and subsequent analysis generates two files, `Exp20_EM350_1_MMStack_Pos0.ome_Left_lvPALM.mat` (Binder), and `Exp20_EM350_1_MMStack_Pos0.ome_RightCali_lvPALM.mat` (tagSrc).
    - They contain the variable `smLinked`, which contains a $M\times 17$ matrix where the first two columns correspond to $x$ and $y$ coordinates in px ($1 \text{ px} = 0.1067\text{ } \mu m$).
    - The `roi.roi` file is generated from drawing an ROI around the cell in [FIJI](https://fiji.sc)
    - These three files need to be paired in the same folder, with a separate folder per cell.
    - For correlation with adhesion images, the file will contain the variables `corrDataStruct` and `focalAdhesionAll`, which additionally have information on tracks in relation to adhesions, and on the adhesions themselves.

2. The diffusional dynamics of these tracks is analyzed using the code in __`diffusional_dynamics`__ (up one level), which relies on [variational Bayes single particle tracking](http://vbspt.sourceforge.net) (vbSPT).
    - In principle, any method to annotate tracks frame-by-frame with slow/fast diffusional states would work.

3. The folders containing the two `.mat` files and `roi.roi` files per each cell are uploaded to Longleaf, as is a folder containing `.mat` files containing the vbSPT-annotated diffusional states. As an example, the file structure might look like, for 3 cells,

        sourcedata/trackdata/C1/Exp20_EM350_1_MMStack_Pos0.ome_Left_lvPALM.mat
        sourcedata/trackdata/C1/Exp20_EM350_1_MMStack_Pos0.ome_RightCali_lvPALM.mat
        sourcedata/trackdata/C1/roi.roi
        sourcedata/trackdata/C2/Exp20_EM350_1_MMStack_Pos0.ome_Left_lvPALM.mat
        sourcedata/trackdata/C2/Exp20_EM350_1_MMStack_Pos0.ome_RightCali_lvPALM.mat
        sourcedata/trackdata/C2/roi.roi
        sourcedata/trackdata/C3/Exp20_EM350_1_MMStack_Pos0.ome_Left_lvPALM.mat
        sourcedata/trackdata/C3/Exp20_EM350_1_MMStack_Pos0.ome_RightCali_lvPALM.mat
        sourcedata/trackdata/C3/roi.roi
        sourcedata/vbspt/source/tag_01.mat
        sourcedata/vbspt/source/tag_02.mat
        sourcedata/vbspt/source/tag_03.mat
        sourcedata/vbspt/source/bin_01.mat
        sourcedata/vbspt/source/bin_02.mat
        sourcedata/vbspt/source/bin_03.mat
        sourcedata/vbspt/source/tag_01.mat
        sourcedata/vbspt/source/tag_02.mat
        sourcedata/vbspt/source/tag_03.mat
        sourcedata/vbspt/source/bin_01.mat
        sourcedata/vbspt/source/bin_02.mat
        sourcedata/vbspt/source/bin_03.mat
        sourcedata/vbspt/result/tag_01_HMManalysis_hidden2.mat
        sourcedata/vbspt/result/tag_02_HMManalysis_hidden2.mat
        sourcedata/vbspt/result/tag_03_HMManalysis_hidden2.mat
        sourcedata/vbspt/result/bin_01_HMManalysis_hidden2.mat
        sourcedata/vbspt/result/bin_02_HMManalysis_hidden2.mat
        sourcedata/vbspt/result/bin_03_HMManalysis_hidden2.mat

4. Either the folder `lvPALM_subroutines` (for non-adhesion-annotated data) or `CDS_subroutines` (for adhesion-annotated date) is uploaded to Longleaf, along with the folder `dependencies`. The folder structure would now look like:
        lvPALM_subroutines/
        dependencies/
        sourcedata/

5. Function(s) to run `cluster_segmentation.m`, and a corresponding Bash file to run them, are uploaded:

*This function calls run cluster_segmentation():*

        function run_analysis()
        % Note -- we keep the dateprefix from the old run (and we might as well keep
        % the parentFigDir) since we need to use the previous refinementII files

        dataloc = 'sourcedata/';
        vbSPTDirLoc = 'sourcedata/vbspt/results/';
        vbSPTMetaDirLoc = 'sourcedata/vbspt/vbspt_source/';
        dateprefix = '20180813';
        parentFigDir = 'figures_20180813/';
        nMinTracks = 10;
        maxRIcutoff = 4;
        numCPUs = 16;
        analysis_mode = 'CDS'; % CDS or lvPALM

        cluster_segmentation(dataloc,dateprefix,parentFigDir,nMinTracks,...
                             maxRIcutoff,vbSPTDirLoc,vbSPTMetaDirLoc,...
                             'temp_01_LLpardir',numCPUs,analysis_mode);`
        end

*This Bash file submits a call to run_analysis() to the SLURM job scheduler on Longleaf*

        #!/bin/bash

        sbatch -p general -N 1 -n 16 -t 60:00:00 --wrap="matlab -nodisplay -nosplash -singleCompThread -r run_analysis\(\) -logfile r01.log"

6. Executing the Bash file should now run the entire pipeline.
7. After this initial analysis, later investigations let us to test several other measurements and apply additional filters. Details on this are provided in the `follow_up_analyses` folder.
