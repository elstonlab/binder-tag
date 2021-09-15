function visualize_f_distribution()
load('annotated_propdists_wRecToClus_minTrk10');

Ti.tag=cell(61,1);
Ti.bin=cell(61,1);
Ti2.tag=cell(61,1);
Ti2.bin=cell(61,1);

for i=1:61
    curr_clus_to_use_tag = ~allanal.propdists.tag{i}.CLVgt2(:) & ...
                            allanal.propdists.tag{i}.startAfter4Sec(:) & ...
                           ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ...
                           ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:);
    curr_clus_to_use_bin = ~allanal.propdists.bin{i}.CLVgt2(:) & ...
                            allanal.propdists.bin{i}.startAfter4Sec(:) & ...
                           ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ...
                           ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:);
                           
    % Try searching only for the 1->0 transition (or 2->0)
    tag_clusters_to_check = find(curr_clus_to_use_tag);
    %Ti.tag{i} = zeros(numel(tag_clusters_to_check),1);
    for j=1:numel(tag_clusters_to_check)
        currclus = tag_clusters_to_check(j);
        transition_type1 = find(diff(logical(allanal.propdists.tag{i}.clus_lifetime_vector{currclus}))==-1);
        Ti.tag{i} = [Ti.tag{i}; diff(transition_type1)*0.02];
    end

    % Try searching only for the 1->0 transition (or 2->0)
    bin_clusters_to_check = find(curr_clus_to_use_bin);
    %Ti.bin{i} = zeros(numel(bin_clusters_to_check),1);
    for j=1:numel(bin_clusters_to_check)
        currclus = bin_clusters_to_check(j);
        transition_type1 = find(diff(logical(allanal.propdists.bin{i}.clus_lifetime_vector{currclus}))==-1);
        Ti.bin{i} = [Ti.bin{i}; diff(transition_type1)*0.02];
    end
    % Try searching only for the 1->0 transition AND the 2->1 transition
    tag_clusters_to_check = find(curr_clus_to_use_tag);
    %Ti2.tag{i} = zeros(numel(tag_clusters_to_check),1);
    for j=1:numel(tag_clusters_to_check)
        currclus = tag_clusters_to_check(j);
        transition_type2 = find(diff((allanal.propdists.tag{i}.clus_lifetime_vector{currclus}))==-1);
        Ti2.tag{i} = [Ti2.tag{i}; diff(transition_type2)*0.02];
    end

    % Try searching only for the 1->0 transition (or 2->0)
    bin_clusters_to_check = find(curr_clus_to_use_bin);
    %Ti2.bin{i} = zeros(numel(bin_clusters_to_check),1);
    for j=1:numel(bin_clusters_to_check)
        currclus = bin_clusters_to_check(j);
        transition_type2 = find(diff((allanal.propdists.bin{i}.clus_lifetime_vector{currclus}))==-1);
        Ti2.bin{i} = [Ti2.bin{i}; diff(transition_type2)*0.02];
    end

end

all_Ti_tag = cell2mat(Ti.tag);
all_Ti_bin = cell2mat(Ti.bin);
all_Ti2_tag = cell2mat(Ti2.tag);
all_Ti2_bin = cell2mat(Ti2.bin);

save('Ti_Ti2','all_Ti_tag','all_Ti_bin','all_Ti2_tag','all_Ti2_bin');

figure('position',[745 151 787 492]);
subplot(2,1,1); hold on;
histogram(cell2mat(Ti.tag),'binedges',0:0.02:5,...
          'normalization','probability');
tag_median = median(cell2mat(Ti.tag));
cylim=get(gca,'ylim');
plot(tag_median*[1,1], cylim, 'k--');
legend('tagSrc clusters',sprintf('median Ti = %.2f', tag_median),...
       'location','northeast');
ylabel('fraction')
subplot(2,1,2); hold on
histogram(cell2mat(Ti.bin),'binedges',0:0.02:5,...
          'normalization','probability',...
          'facecolor','r');
bin_median = median(cell2mat(Ti.bin));
cylim=get(gca,'ylim');
plot(bin_median*[1,1], cylim, 'k--');
xlabel('Ti')
ylabel('fraction')
legend('binder clusters',sprintf('median Ti = %.2f', bin_median),...
       'location','northeast');
savefig(['Ti_distribution.fig']);
print(gcf,'-dpng',['Ti_distribution.png'],'-r300');
print(gcf,'-dsvg','-painters',['Ti_distribution.svg'],'-r300');


figure('position',[745 151 787 492]);
subplot(2,1,1); hold on;
histogram(cell2mat(Ti2.tag),'binedges',0:0.02:5,...
          'normalization','probability');
tag_median = median(cell2mat(Ti2.tag));
cylim=get(gca,'ylim');
plot(tag_median*[1,1], cylim, 'k--');
legend('tagSrc clusters',sprintf('median Ti2 = %.2f', tag_median),...
       'location','northeast');
ylabel('fraction')
subplot(2,1,2); hold on
histogram(cell2mat(Ti2.bin),'binedges',0:0.02:5,...
          'normalization','probability',...
          'facecolor','r');
bin_median = median(cell2mat(Ti2.bin));
cylim=get(gca,'ylim');
plot(bin_median*[1,1], cylim, 'k--');
xlabel('Ti2')
ylabel('fraction')
legend('binder clusters',sprintf('median Ti2 = %.2f', bin_median),...
       'location','northeast');
savefig(['Ti2_distribution.fig']);
print(gcf,'-dpng',['Ti2_distribution.png'],'-r300');
print(gcf,'-dsvg','-painters',['Ti2_distribution.svg'],'-r300');
end