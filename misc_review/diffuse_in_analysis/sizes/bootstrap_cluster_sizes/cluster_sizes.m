function cluster_sizes()
load('C:\Users\mike\Documents\GitHub\binder-tag\initial_review\diffuse_in_analysis\annotated_propdists_wRecToClus_minTrk10.mat')
ar50_tag = cell2mat(cellfun(@(x) x.allLoc.KDE.area50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
ar50_bin = cell2mat(cellfun(@(x) x.allLoc.KDE.area50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));

[ar50_tag_bs,tmp]=bootstrp(100000,@nanmedian,ar50_tag);
[ar50_bin_bs,tmp]=bootstrp(100000,@nanmedian,ar50_bin);

data = [ar50_tag_bs,ar50_bin_bs];

writematrix(data,'tagsrc_binder_ar50_removeDiffuseIn_bootstrapped.csv');

end
