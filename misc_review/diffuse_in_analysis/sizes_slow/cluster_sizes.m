load('C:\Users\mike\Documents\GitHub\binder-tag\initial_review\diffuse_in_analysisannotated_propdists_wRecToClus.mat')

slowar50_tag = cell2mat(cellfun(@(x) x.vbSPT.slowarea50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
slowar50_bin = cell2mat(cellfun(@(x) x.vbSPT.slowarea50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));
fastar50_tag = cell2mat(cellfun(@(x) x.vbSPT.fastarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
fastar50_bin = cell2mat(cellfun(@(x) x.vbSPT.fastarea50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));

figure('position',[705   912   303   216])
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
