function twotrackssimultaneous()
% Visualize proprotion of clusters with >2 tracks simultaneously.
% Within that proportion, what proportion of the cluster's lifetime has two or more?
load('annotated_propdists_wRecToClus_minTrk10');

is_gt2.tag=cell(61,1);
is_gt2.bin=cell(61,1);
time_gt2.tag=cell(61,1);
time_gt2.bin=cell(61,1);

for i=1:61
    curr_clus_to_use_tag = ~allanal.propdists.tag{i}.CLVgt2(:) & ...
                            allanal.propdists.tag{i}.startAfter4Sec(:) & ...
                           ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ...
                           ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:);
    curr_clus_to_use_bin = ~allanal.propdists.bin{i}.CLVgt2(:) & ...
                            allanal.propdists.bin{i}.startAfter4Sec(:) & ...
                           ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ...
                           ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:);
    is_gt2.tag{i} = cellfun(@(x) any(x>=2),...
                        allanal.propdists.tag{i}.clus_lifetime_vector(curr_clus_to_use_tag));
    is_gt2.bin{i} = cellfun(@(x) any(x>=2),...
                        allanal.propdists.bin{i}.clus_lifetime_vector(curr_clus_to_use_bin));
    tmptag = allanal.propdists.tag{i}.clus_lifetime_vector(curr_clus_to_use_tag);
    tmpbin = allanal.propdists.bin{i}.clus_lifetime_vector(curr_clus_to_use_bin);
                    
    time_gt2.tag{i} = cellfun(@(x) sum(x>=2)/numel(x),...
                              tmptag(is_gt2.tag{i}));
    time_gt2.bin{i} = cellfun(@(x) sum(x>=2)/numel(x),...
                              tmpbin(is_gt2.bin{i}));

                    
end
is_gt2_mat_tag = cell2mat(is_gt2.tag);
is_gt2_mat_bin = cell2mat(is_gt2.bin);


time_gt2_mat_tag = cell2mat(time_gt2.tag);
time_gt2_mat_bin = cell2mat(time_gt2.bin);

sum(is_gt2_mat_tag)/numel(is_gt2_mat_tag)*100
sum(is_gt2_mat_bin)/numel(is_gt2_mat_bin)*100

figure('position',[484 318 273 492]);
subplot(2,1,1); hold on;
histogram(time_gt2_mat_tag,'binedges',0:0.05:1,...
          'normalization','probability');
tag_median = median(time_gt2_mat_tag);
cylim=get(gca,'ylim');
plot(tag_median*[1,1], cylim, 'k--');
legend('tagSrc clusters',sprintf('median f = %.2f', tag_median),...
       'location','northeast');
ylabel('fraction')
subplot(2,1,2); hold on
histogram(time_gt2_mat_bin,'binedges',0:0.05:1,...
          'normalization','probability',...
          'facecolor','r');
bin_median = median(time_gt2_mat_bin);
cylim=get(gca,'ylim');
plot(bin_median*[1,1], cylim, 'k--');
xlabel('fraction lifetime spent with 2 tracks present')
ylabel('fraction')
legend('binder clusters',sprintf('median f = %.2f', bin_median),...
       'location','northeast');
savefig(['fraction_lifetime_2tracks.fig']);
print(gcf,'-dpng',['fraction_lifetime_2tracks.png'],'-r300');
print(gcf,'-dsvg','-painters',['fraction_lifetime_2tracks.svg'],'-r300');
end