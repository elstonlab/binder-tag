load('I:\Elston Lab\MikePablo\20180504-bindertag\20190617_remeasuring_properties\lifetimes\track_lifetime_analysis.mat')
lifetimes_tag_clus = cell2mat(tag_cluslifetimes)*.02;
lifetimes_tag_nonclus = cell2mat(tag_noncluslifetimes)*.02;
lifetimes_bin_clus = cell2mat(bin_cluslifetimes)*.02;
lifetimes_bin_nonclus = cell2mat(bin_noncluslifetimes)*.02;

[LT_tag_clus_bs,~]=bootstrp(100,@median,lifetimes_tag_clus);
[LT_tag_nonclus_bs,~]=bootstrp(100,@median,lifetimes_tag_nonclus);
[LT_bin_clus_bs,~]=bootstrp(100,@median,lifetimes_bin_clus);
[LT_bin_nonclus_bs,~]=bootstrp(100,@median,lifetimes_bin_nonclus);

figure;
subplot(1,2,1); hold on;
histogram(LT_tag_clus_bs,'binedges',(0:.0200:2),'facecolor',[51 102 51]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(LT_tag_nonclus_bs,'binedges',(0:.0200:2),'edgecolor','k','displaystyle','stairs','normalization','pdf');
xlabel('track lifetime (s)')
ylabel('fraction of tagSrc tracks');
legend('clustered','non-clustered');
set(gca,'fontsize',14);

subplot(1,2,2); hold on;
histogram(LT_bin_clus_bs,'binedges',(0:.0200:2),'facecolor',[255 51 51]/255,'edgecolor','none','normalization','pdf');
ylabel('probability density')
histogram(LT_bin_nonclus_bs,'binedges',(0:.0200:2),'edgecolor','k','displaystyle','stairs','normalization','pdf');
xlabel('track lifetime (s)')
ylabel('fraction of tagSrc tracks');
legend('clustered','non-clustered');
set(gca,'fontsize',14);