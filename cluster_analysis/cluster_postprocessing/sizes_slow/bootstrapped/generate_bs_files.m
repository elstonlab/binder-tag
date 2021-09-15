load('J:\Elston Lab\MikePablo\20180504-bindertag\20190528_rectoclus\annotated_propdists_wRecToClus.mat')
slowar50_tag = cell2mat(cellfun(@(x) x.vbSPT.slowarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
slowar50_bin = cell2mat(cellfun(@(x) x.vbSPT.slowarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));
fastar50_tag = cell2mat(cellfun(@(x) x.vbSPT.fastarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
fastar50_bin = cell2mat(cellfun(@(x) x.vbSPT.fastarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));


% add bootstrapping
[slowar50_tag_bs,~]=bootstrp(100000,@nanmedian,slowar50_tag);
[slowar50_bin_bs,~]=bootstrp(100000,@nanmedian,slowar50_bin);
[fastar50_tag_bs,~]=bootstrp(100000,@nanmedian,fastar50_tag);
[fastar50_bin_bs,~]=bootstrp(100000,@nanmedian,fastar50_bin);

data = [slowar50_tag_bs,slowar50_bin_bs,...
        fastar50_tag_bs,fastar50_bin_bs];

writematrix(data,'slowTS-slowB-fastTS-fastB_ar50_bootstrapped.csv');

paired_diff50_tag = fastar50_tag - slowar50_tag;
paired_diff50_bin = fastar50_bin - slowar50_bin;

figure('position',[204 725 554 222]);
hold on;
histogram(paired_diff50_tag,'binedges',-2e5:0.05e5:2e5,'facecolor',[131 161 131]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(paired_diff50_bin,'binedges',-2e5:0.05e5:2e5,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('area, paired differences (fast-slow) (nm2)')
savefig('paired_difference');
print('-dsvg','-painters','paired_difference.svg','-r300');
print('-dpng','paired_difference.png','-r300');