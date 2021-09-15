load('J:\Elston Lab\MikePablo\20180504-bindertag\20190528_rectoclus\annotated_propdists_wRecToClus.mat')
slowar50_tag = cell2mat(cellfun(@(x) x.vbSPT.slowarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
slowar50_bin = cell2mat(cellfun(@(x) x.vbSPT.slowarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));
fastar50_tag = cell2mat(cellfun(@(x) x.vbSPT.fastarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
fastar50_bin = cell2mat(cellfun(@(x) x.vbSPT.fastarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));

figure;
hold on;
histogram(slowar50_tag,'binedges',0:0.05e5:2.5e5,'edgecolor','k','displaystyle','stairs','normalization','pdf','linestyle','--');
ylabel('probability density')
histogram(slowar50_bin,'binedges',0:0.05e5:2.5e5,'edgecolor','r','displaystyle','stairs','normalization','pdf','linestyle','--');
xlabel('slow area (nm2)')
set(gca,'fontsize',14,'ylim',[0 3e-5],'ytick',[0 1e-5 2e-5 3e-5]);

histogram(fastar50_tag,'binedges',0:0.05e5:2.5e5,'edgecolor','k','displaystyle','stairs','normalization','pdf');
ylabel('probability density')
histogram(fastar50_bin,'binedges',0:0.05e5:2.5e5,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('fast area (nm2)')
set(gca,'fontsize',14,'ylim',[0 3e-5],'ytick',[0 1e-5 2e-5 3e-5]);
legend('tagSrc, slow','Binder, slow','tagSrc, fast','Binder, fast');