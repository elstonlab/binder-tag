load('nMinTrack_result_10/annotated_propdists_wRecToClus_minTrk10.mat')

tag_CLVgt2_percent = zeros(61,1);
bin_CLVgt2_percent = zeros(61,1);

for i=1:61
    boolidx_tag = allanal.propdists.tag{i}.startAfter4Sec(:) & ...
                  ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ...
                  ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:);
    boolidx_bin = allanal.propdists.bin{i}.startAfter4Sec(:) & ...
                  ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ...
                  ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:);
    tag_CLVgt2_percent(i) = 100*sum(allanal.propdists.tag{i}.CLVgt2(boolidx_tag)) /  ...
                              numel(allanal.propdists.tag{i}.CLVgt2(boolidx_tag));
    bin_CLVgt2_percent(i) = 100*sum(allanal.propdists.bin{i}.CLVgt2(boolidx_bin)) /  ...
                              numel(allanal.propdists.bin{i}.CLVgt2(boolidx_bin));
end
figure;
histogram(tag_CLVgt2_percent,'binedges',0:10:100)
xlabel({'Percentage of tagSrc clusters with','>2 tracks present at some point'})
ylabel('# cells')
figure;
histogram(bin_CLVgt2_percent,'binedges',0:10:100)
xlabel({'Percentage of Binder clusters with','>2 tracks present at some point'})
ylabel('# cells')
