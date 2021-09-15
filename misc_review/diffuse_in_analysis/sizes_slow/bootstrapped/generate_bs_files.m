function generate_bs_files()
load('C:\Users\mike\Documents\GitHub\binder-tag\initial_review\diffuse_in_analysis\annotated_propdists_wRecToClus_minTrk10.mat')
slowar50_tag = cell2mat(cellfun(@(x) x.vbSPT.slowarea50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
slowar50_bin = cell2mat(cellfun(@(x) x.vbSPT.slowarea50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));
fastar50_tag = cell2mat(cellfun(@(x) x.vbSPT.fastarea50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
fastar50_bin = cell2mat(cellfun(@(x) x.vbSPT.fastarea50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));


% add bootstrapping
[slowar50_tag_bs,~]=bootstrp(100000,@nanmedian,slowar50_tag);
[slowar50_bin_bs,~]=bootstrp(100000,@nanmedian,slowar50_bin);
[fastar50_tag_bs,~]=bootstrp(100000,@nanmedian,fastar50_tag);
[fastar50_bin_bs,~]=bootstrp(100000,@nanmedian,fastar50_bin);

data = [slowar50_tag_bs,slowar50_bin_bs,...
        fastar50_tag_bs,fastar50_bin_bs];

writematrix(data,'slowTS-slowB-fastTS-fastB_ar50_bootstrapped_removeDiffuseIn.csv');

paired_diff50_tag = fastar50_tag - slowar50_tag;
paired_diff50_bin = fastar50_bin - slowar50_bin;

figure('position',[349   419   303   244])
subplot(1,2,1); hold on;
histogram(slowar50_tag,'binedges',0:0.05e5:2.5e5,'facecolor',[51 204 255]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(slowar50_bin,'binedges',0:0.05e5:2.5e5,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('slow area (nm2)')
set(gca,'fontsize',14,'ylim',[0 3e-5],'ytick',[0 1e-5 2e-5 3e-5]);
text(1.75e5,2.7e-5,sprintf('n(tagSrc) = %i',numel(slowar50_tag)),'horizontalAlignment','center','color',[51 204 255]/255,'fontsize',10);
text(1.75e5,2.2e-5,sprintf('n(Binder) = %i',numel(slowar50_bin)),'horizontalAlignment','center','color','r','fontsize',10);
subplot(1,2,2); hold on;
histogram(fastar50_tag,'binedges',0:0.05e5:2.5e5,'facecolor',[51 204 255]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(fastar50_bin,'binedges',0:0.05e5:2.5e5,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('fast area (nm2)')
set(gca,'fontsize',14,'ylim',[0 3e-5],'ytick',[0 1e-5 2e-5 3e-5]);
text(1.75e5,2.7e-5,sprintf('n(tagSrc) = %i',numel(fastar50_tag)),'horizontalAlignment','center','color',[51 204 255]/255,'fontsize',10);
text(1.75e5,2.2e-5,sprintf('n(Binder) = %i',numel(fastar50_bin)),'horizontalAlignment','center','color','r','fontsize',10);
savefig('slowfast_removeDiffuseIn');
print('-dsvg','-painters','slowfast_removeDiffuseIn.svg','-r300');
print('-dpng','slowfast_removeDiffuseIn.png','-r300');

figure('position',[349   419   303   244])
hold on;
histogram(paired_diff50_tag,'binedges',-2e5:0.05e5:2e5,'facecolor',[131 161 131]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(paired_diff50_bin,'binedges',-2e5:0.05e5:2e5,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('area, paired differences (fast-slow) (nm2)')
savefig('paired_difference_removeDiffuseIn');
print('-dsvg','-painters','paired_difference_removeDiffuseIn.svg','-r300');
print('-dpng','paired_difference_removeDiffuseIn.png','-r300');
end