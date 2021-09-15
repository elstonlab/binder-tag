function stack_clusterNtrack_distributions_4_to_9()
ntracks=4:9;
figure('position',[127   683   653   200]);
summarypath = 'C:\\Users\\mike\\Documents\\GitHub\\binder-tag\\initial_review\\vary_nMinTrk\\cluster_lifetime_and_Ntracks_summary.mat';
for i=1:numel(ntracks)
    curr_ntrack = ntracks(i);
    sourcepath = sprintf('C:\\Users\\mike\\Documents\\GitHub\\binder-tag\\initial_review\\vary_nMinTrk\\nMinTrack_result_%02d\\annotated_propdists_wRecToClus_minTrk%02d',...
                         curr_ntrack, curr_ntrack);
    outpath = sprintf('C:\\Users\\mike\\Documents\\GitHub\\binder-tag\\initial_review\\vary_nMinTrk\\nMinTrack_result_%02d\\', curr_ntrack);
    
    visualize_datasets(sourcepath, outpath);
    
    
end
lgd=legend('4','5','6','7','8','9');
lgd.FontSize=8;
savefig(['stack_nTrack_distributions_4_to_9.fig']);
print(gcf,'-dpng',['stack_nTrack_distributions_4_to_9.png'],'-r300');
print(gcf,'-dsvg','-painters',['stack_nTrack_distributions_4_to_9.svg'],'-r300');

end

function [LT_tag, LT_bin] = visualize_datasets(sourcepath, outpath)
load(sourcepath)
N_tag = cell2mat(cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false));
N_bin = cell2mat(cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false));
[f_tag,xi_tag]=ksdensity(N_tag, 0:1:60);
[f_bin,xi_bin]=ksdensity(N_bin, 0:1:60);

subplot(1,2,1); hold on
ylabel('pdf');
plot(xi_tag, f_tag,'linewidth',2);
xlabel('tagSrc tracks per cluster')
set(gca,'fontsize',14,'ylim',[0 0.20],'ytick',0:.05:.20,'xlim',[0 50]);
subplot(1,2,2); hold on
plot(xi_bin, f_bin,'linewidth',2);
xlabel('binder tracks per cluster')
set(gca,'fontsize',14,'ylim',[0 0.20],'ytick',0:.05:.20,'xlim',[0 50]);
end
