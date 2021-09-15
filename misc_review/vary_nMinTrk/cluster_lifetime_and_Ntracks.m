function cluster_lifetime_and_Ntracks()
ntracks=3:20;
summarypath = 'C:\\Users\\mike\\Documents\\GitHub\\binder-tag\\initial_review\\vary_nMinTrk\\cluster_lifetime_and_Ntracks_summary.mat';
for i=1:numel(ntracks)
    curr_ntrack = ntracks(i);
    sourcepath = sprintf('C:\\Users\\mike\\Documents\\GitHub\\binder-tag\\initial_review\\vary_nMinTrk\\nMinTrack_result_%02d\\annotated_propdists_wRecToClus_minTrk%02d',...
                         curr_ntrack, curr_ntrack);
    outpath = sprintf('C:\\Users\\mike\\Documents\\GitHub\\binder-tag\\initial_review\\vary_nMinTrk\\nMinTrack_result_%02d\\', curr_ntrack);
    
    [LT_tag_median(i), LT_bin_median(i), N_tag_median(i), N_bin_median(i),...
     lifetime_tag{i}, lifetime_bin{i}, N_tag{i}, N_bin{i}] = process_dataset(sourcepath, outpath);
    
    
end
save(summarypath, 'LT_tag_median','LT_bin_median','N_tag_median','N_bin_median',...
                   'lifetime_tag', 'lifetime_bin', 'N_tag', 'N_bin', 'ntracks');
end

function [LT_tag_median, LT_bin_median, N_tag_median, N_bin_median,...
          lifetime_tag, lifetime_bin, N_tag, N_bin] = process_dataset(sourcepath, outpath)
load(sourcepath)

lifetime_tag = cell2mat(cellfun(@(x) x.lifetime(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
lifetime_bin = cell2mat(cellfun(@(x) x.lifetime(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));

N_tag = cell2mat(cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
N_bin = cell2mat(cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));

figure('position',[187   379   303   216]); hold on
histogram(lifetime_tag,'binedges',0:1:60,'facecolor',[51 204 255]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(lifetime_bin,'binedges',0:1:60,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('lifetime (s)')
set(gca,'fontsize',14,'ylim',[0 0.20],'ytick',0:.05:.20);

savefig([outpath, 'cluster_lifetime.fig']);
print(gcf,'-dpng',[outpath, 'cluster_lifetime.png'],'-r300');

figure('position',[187   379   303   216]); hold on
histogram(N_tag,'binedges',0:1:60,'facecolor',[51 204 255]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(N_bin,'binedges',0:1:60,'edgecolor','r','displaystyle','stairs','normalization','pdf');
xlabel('Tracks per cluster')
set(gca,'fontsize',14,'ylim',[0 0.20],'ytick',0:.05:.20);
savefig([outpath, 'cluster_Ntracks.fig']);
print(gcf,'-dpng',[outpath, 'cluster_Ntracks.png'],'-r300');

LT_tag_median = nanmedian(lifetime_tag);
LT_bin_median = nanmedian(lifetime_bin);
N_tag_median = nanmedian(N_tag);
N_bin_median = nanmedian(N_bin);
end
