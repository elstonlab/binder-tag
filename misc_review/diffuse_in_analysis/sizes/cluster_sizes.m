load('C:\Users\mike\Documents\GitHub\binder-tag\initial_review\diffuse_in_analysis\annotated_propdists_wRecToClus_minTrk10.mat')
ar50_tag = cell2mat(cellfun(@(x) x.allLoc.KDE.area50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
ar50_bin = cell2mat(cellfun(@(x) x.allLoc.KDE.area50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));

figure('position',[349   419   303   244])
hold on;
histogram(ar50_tag,'binedges',0:0.05e5:2.5e5,'facecolor',[51 204 255]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(ar50_bin,'binedges',0:0.05e5:2.5e5,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('area (nm2)')
set(gca,'fontsize',14,'ylim',[0 3e-5],'ytick',[0 1e-5 2e-5 3e-5]);
text(1.75e5,2.7e-5,sprintf('n(tagSrc) = %i',numel(ar50_tag)),'horizontalAlignment','center','color',[51 204 255]/255,'fontsize',10);
text(1.75e5,2.2e-5,sprintf('n(Binder) = %i',numel(ar50_bin)),'horizontalAlignment','center','color','r','fontsize',10);

% Compare vs removed diffusing in
ar50_tag = cell2mat(cellfun(@(x) x.allLoc.KDE.area50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
ar50_bin = cell2mat(cellfun(@(x) x.allLoc.KDE.area50_removeDiffuseIn(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));

figure('position',[349   419   303   244])
hold on;
histogram(ar50_tag,'binedges',0:0.05e5:2.5e5,'facecolor',[51 204 255]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(ar50_bin,'binedges',0:0.05e5:2.5e5,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('area, no diffusing in (nm2)')
set(gca,'fontsize',14,'ylim',[0 7e-5],'ytick',[0:1e-5:7e-5]);
text(1.75e5,5.7e-5,sprintf('n(tagSrc) = %i',numel(ar50_tag)),'horizontalAlignment','center','color',[51 204 255]/255,'fontsize',10);
text(1.75e5,5.2e-5,sprintf('n(Binder) = %i',numel(ar50_bin)),'horizontalAlignment','center','color','r','fontsize',10);
