function median_lifetime_and_RI_by_ntracks()
load('annotated_propdists_wRecToClus_minTrk10.mat')

track_groupings = 10:30;

n_track_groupings = numel(track_groupings);

rec_int_by_cluster.tag=cell(61,n_track_groupings);
rec_int_by_cluster.bin=cell(61,n_track_groupings);
lifetime.tag=cell(61,n_track_groupings);
lifetime.bin=cell(61,n_track_groupings);

for j=1:n_track_groupings
    for i=1:61
        
        curr_clus_to_use_tag = ~allanal.propdists.tag{i}.CLVgt2(:) & ...
                                allanal.propdists.tag{i}.startAfter4Sec(:) & ...
                               ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ...
                               ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:) & ...
                                allanal.propdists.tag{i}.N == track_groupings(j);
        curr_clus_to_use_bin = ~allanal.propdists.bin{i}.CLVgt2(:) & ...
                                allanal.propdists.bin{i}.startAfter4Sec(:) & ...
                               ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ...
                               ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:) & ...
                                allanal.propdists.bin{i}.N == track_groupings(j);

        
        rec_int_by_cluster.tag{i,j} = cell2mat(allanal.propdists.tag{i}.rec_int_by_cluster(curr_clus_to_use_tag));
        rec_int_by_cluster.bin{i,j} = cell2mat(allanal.propdists.bin{i}.rec_int_by_cluster(curr_clus_to_use_bin));
        lifetime.tag{i,j} = allanal.propdists.tag{i}.lifetime(curr_clus_to_use_tag);
        lifetime.bin{i,j} = allanal.propdists.bin{i}.lifetime(curr_clus_to_use_bin);


    end
end
RI_median_tag = zeros(n_track_groupings,1);
RI_median_bin = zeros(n_track_groupings,1);
lifetime_median_tag = zeros(n_track_groupings,1);
lifetime_median_bin = zeros(n_track_groupings,1);
lifetime_ksdensity_tag = cell(n_track_groupings,1);
lifetime_ksdensity_bin = cell(n_track_groupings,1);
for i=1:n_track_groupings
    RI_median_tag(i) = median(cell2mat(rec_int_by_cluster.tag(:,i)));
    RI_median_bin(i) = median(cell2mat(rec_int_by_cluster.bin(:,i)));
    
    lifetime_median_tag(i) = median(cell2mat(lifetime.tag(:,i)));
    lifetime_median_bin(i) = median(cell2mat(lifetime.bin(:,i)));
    
    lifetime_ksdensity_tag{i} = ksdensity(cell2mat(lifetime.tag(:,i)),0:.1:60);
    lifetime_ksdensity_bin{i} = ksdensity(cell2mat(lifetime.bin(:,i)),0:.1:60);
end

figure('position',[2331 510 389 296]); 
hold on;
plot(track_groupings,lifetime_median_tag,'.')
plot(track_groupings,lifetime_median_bin,'r.')
xlabel('# tracks in cluster')
ylabel('lifetime median (s)')
ylabel('cluster lifetime median (s)')
legend('tagSrc','Binder','location','northwest')
savefig(['by_track_medianClusLifetime.fig']);
print(gcf,'-dpng',['by_track_medianClusLifetime.png'],'-r300');
print(gcf,'-dsvg','-painters',['by_track_medianClusLifetime.svg'],'-r300');


figure('position',[2331 510 389 296]); 
hold on;
plot(track_groupings,RI_median_tag,'.')
plot(track_groupings,RI_median_bin,'r.')
xlabel('# tracks in cluster')
ylabel('lifetime median (s)')
ylabel('recruitment interval median (s)')
legend('tagSrc','Binder','location','northeast')
savefig(['by_track_medianRI.fig']);
print(gcf,'-dpng',['by_track_medianRI.png'],'-r300');
print(gcf,'-dsvg','-painters',['by_track_medianRI.svg'],'-r300');

cmap=copper(n_track_groupings);
figure('position',[857   441   364   405]); 
subplot(2,1,1)
hold on
for i=1:20
    plot(0:.1:60,lifetime_ksdensity_tag{i},'color',cmap(i,:));
end
ylabel('density')
subplot(2,1,2)
hold on
for i=1:n_track_groupings
    plot(0:.1:60,lifetime_ksdensity_bin{i},'color',cmap(i,:));
end
ylabel('density')
xlabel('cluster lifetime')
subplot(2,1,1)
toplgd=legend(arrayfun(@num2str, track_groupings, 'UniformOutput', false),'location','eastoutside');
subplot(2,1,2)
botlgd=legend(arrayfun(@num2str, track_groupings, 'UniformOutput', false),'location','eastoutside');
botlgd.Visible='off';
savefig(['by_track_RI_distributions.fig']);
print(gcf,'-dpng',['by_track_RI_distributions.png'],'-r300');
print(gcf,'-dsvg','-painters',['by_track_RI_distributions.svg'],'-r300');
end