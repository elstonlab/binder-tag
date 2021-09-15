function visualize_f_distribution()
load('annotated_propdists_wRecToClus_minTrk10');

f.tag=cell(61,1);
f.bin=cell(61,1);

for i=1:61
    curr_clus_to_use_tag = ~allanal.propdists.tag{i}.CLVgt2(:) & ...
                            allanal.propdists.tag{i}.startAfter4Sec(:) & ...
                           ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ...
                           ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:);
    curr_clus_to_use_bin = ~allanal.propdists.bin{i}.CLVgt2(:) & ...
                            allanal.propdists.bin{i}.startAfter4Sec(:) & ...
                           ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ...
                           ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:);
                           
    f.tag{i} = cellfun(@(x) sum(x>=1)/numel(x),...
                       allanal.propdists.tag{i}.clus_lifetime_vector(curr_clus_to_use_tag));
    f.bin{i} = cellfun(@(x) sum(x>=1)/numel(x),...
                       allanal.propdists.bin{i}.clus_lifetime_vector(curr_clus_to_use_bin));
                    
end
figure('position',[484 318 273 492]);
subplot(2,1,1); hold on;
histogram(cell2mat(f.tag),'binedges',0:0.05:1,...
          'normalization','probability');
tag_median = median(cell2mat(f.tag));
cylim=get(gca,'ylim');
plot(tag_median*[1,1], cylim, 'k--');
legend('tagSrc clusters',sprintf('median f = %.2f', tag_median),...
       'location','northwest');
ylabel('fraction')
subplot(2,1,2); hold on
histogram(cell2mat(f.bin),'binedges',0:0.05:1,...
          'normalization','probability',...
          'facecolor','r');
bin_median = median(cell2mat(f.bin));
cylim=get(gca,'ylim');
plot(bin_median*[1,1], cylim, 'k--');
xlabel('f')
ylabel('fraction')
legend('binder clusters',sprintf('median f = %.2f', bin_median),...
       'location','northwest');
savefig(['f_distribution.fig']);
print(gcf,'-dpng',['f_distribution.png'],'-r300');
print(gcf,'-dsvg','-painters',['f_distribution.svg'],'-r300');
end