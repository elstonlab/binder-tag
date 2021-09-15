function cluster_sizes()
load('J:\Elston Lab\MikePablo\20180504-bindertag\20190528_rectoclus\annotated_propdists_wRecToClus.mat')
ar50_tag = cell2mat(cellfun(@(x) x.allLoc.KDE.area50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
ar50_bin = cell2mat(cellfun(@(x) x.allLoc.KDE.area50(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));
ar95_tag = cell2mat(cellfun(@(x) x.allLoc.KDE.area95(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
ar95_bin = cell2mat(cellfun(@(x) x.allLoc.KDE.area95(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));

lifetime_tag = cell2mat(cellfun(@(x) x.lifetime(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
lifetime_bin = cell2mat(cellfun(@(x) x.lifetime(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));

% [lt_tag_bs,tmp]=bootstrp(100000,@median,lifetime_tag);
% [lt_bin_bs,tmp]=bootstrp(100000,@median,lifetime_bin);

[ar50_tag_bs,~]=bootstrp(100000,@median,ar50_tag);
[ar50_bin_bs,~]=bootstrp(100000,@median,ar50_bin);
% [ar95_tag_bs,tmp]=bootstrp(100000,@median,ar95_tag);
% [ar95_bin_bs,tmp]=bootstrp(100000,@median,ar95_bin);

figure('position',[705   912   303   216])
subplot(1,2,1); hold on;
histogram(ar50_tag,'binedges',0:0.05e5:2.5e5,'facecolor',[51 204 255]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(ar50_bin,'binedges',0:0.05e5:2.5e5,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('area (nm2)')
set(gca,'fontsize',14,'ylim',[0 3e-5],'ytick',[0 1e-5 2e-5 3e-5]);
text(1.75e5,2.7e-5,sprintf('n(tagSrc) = %i',numel(ar50_tag)),'horizontalAlignment','center','color',[51 204 255]/255,'fontsize',10);
text(1.75e5,2.2e-5,sprintf('n(Binder) = %i',numel(ar50_tag)),'horizontalAlignment','center','color','r','fontsize',10);

subplot(1,2,2); hold on
histogram(lifetime_tag,'binedges',0:1:60,'facecolor',[51 204 255]/255,'edgecolor','none','normalization','pdf');
%ylabel('probability density')
histogram(lifetime_bin,'binedges',0:1:60,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('lifetime (s)')
set(gca,'fontsize',14,'ylim',[0 0.15],'ytick',0:.05:.15);
end
