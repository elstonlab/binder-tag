load('C:\Users\mike\Documents\GitHub\binder-tag\initial_review\CAtagSrc_vs_Binder_analysis\20201212\followupanalysis\annotated_propdists_wRecToClus.mat')

curr_filelocation = 'C:\\Users\\mike\\Documents\\GitHub\\binder-tag\\initial_review\\CAtagSrc_vs_Binder_analysis\\20201212\\results\\inputdata\\vbspt_source\\';
N_tag = cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false);
N_bin = cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false);

percent_tracks_in_clus_tag = zeros(15,1);
percent_tracks_in_clus_bin = zeros(15,1);


% For each cell
for i=1:15
    % # tracks in clusters after stringent filtering.
    curr_N_tag = sum(N_tag{i});
    curr_N_bin = sum(N_bin{i});
    
    curr_file_index = i;
    tag_track_IDs = load(sprintf('%stag_%02d.mat',curr_filelocation,curr_file_index),'vbspt_trackIDs');
    bin_track_IDs = load(sprintf('%sbin_%02d.mat',curr_filelocation,curr_file_index),'vbspt_trackIDs');
    
    N_usable_tag = numel(tag_track_IDs.vbspt_trackIDs);
    N_usable_bin = numel(bin_track_IDs.vbspt_trackIDs);
    
    percent_tracks_in_clus_tag(i) = curr_N_tag/N_usable_tag*100;
    percent_tracks_in_clus_bin(i) = curr_N_bin/N_usable_bin*100;    
end
% Get the total # tracks in clusters after filtering -- use the annotated_propdists data


figure('position',[705   912   303   216]); hold on;
plotSpread([percent_tracks_in_clus_tag,percent_tracks_in_clus_bin],'distributionColors',[51 102 51; 255 51 51]/255);
ylabel('tracks in clusters (%)')
set(gca,'fontsize',14,'xtick',[1 2],'xticklabel',{'tagSrc-CA','Binder'},'ylim',[0 40]);
%text(1.75e5,2.7e-5,sprintf('n(tagSrc) = %i',numel(ar50_tag)),'horizontalAlignment','center','color',[51 204 255]/255,'fontsize',10);
%text(1.75e5,2.2e-5,sprintf('n(Binder) = %i',numel(ar50_tag)),'horizontalAlignment','center','color','r','fontsize',10);
print(gcf,'-dpng','tracks_in_clusters.png','-r300');
print(gcf,'-dsvg','-painters','tracks_in_clusters.svg','-r300');