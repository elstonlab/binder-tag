function measure_time_in_slow_state()

if ~exist('track_slowFrac_time.mat')

filelocations = {'C:\\Users\\mike\\Documents\\GitHub\\binder-tag\\initial_review\\CAtagSrc_vs_Binder_analysis\\20201212\\results\\inputdata\\'};

% We need to generate the expdata files to get track lifetimes..
all_path_1.tag = {};
all_path_1.bin = {};
all_file_1.tag = {};
all_file_1.bin = {};
for i=1:1
    [curr_path,curr_file] = batchGetPath(filelocations{1},'RightCali_lvPALM','mat');
    all_path_1.tag = [all_path_1.tag;curr_path];
    all_file_1.tag = [all_file_1.tag;curr_file];
    [curr_path,curr_file] = batchGetPath(filelocations{1},'Left_lvPALM','mat');
    all_path_1.bin = [all_path_1.bin;curr_path];
    all_file_1.bin = [all_file_1.bin;curr_file];
end
assert(numel(all_path_1.tag)==numel(all_path_1.bin));

annotated_propdist=load('C:\Users\mike\Documents\GitHub\binder-tag\initial_review\CAtagSrc_vs_Binder_analysis\20201212\followupanalysis\annotated_propdists_wRecToClus.mat','allanal');

N_tag = cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),annotated_propdist.allanal.propdists.tag,'uniformoutput',false);
N_bin = cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),annotated_propdist.allanal.propdists.bin,'uniformoutput',false);

ncells=15;

percent_tracks_in_clus_tag = zeros(ncells,1);
percent_tracks_in_clus_bin = zeros(ncells,1);


meas1 = load('C:\\Users\\mike\\Documents\\GitHub\\binder-tag\\initial_review\\CAtagSrc_vs_Binder_analysis\\20201212\\results\\20201215measured_cluster_props.mat');
     

tag_clus_totSlowFrames = zeros(ncells,1);
bin_clus_totSlowFrames = zeros(ncells,1);
tag_nonclus_totSlowFrames = zeros(ncells,1);
bin_nonclus_totSlowFrames = zeros(ncells,1);

tag_clus_totFrames = zeros(ncells,1);
bin_clus_totFrames = zeros(ncells,1);
tag_nonclus_totFrames = zeros(ncells,1);
bin_nonclus_totFrames = zeros(ncells,1);
    
% For each cell
for i=1:ncells
    fprintf('on cell %i\n',i);
    % Get the total usable tracks (read the vbSPT source file and check the # of vbspt_trackIDs.

    curr_filelocation = filelocations{1};
    curr_file_index = i;

    curr_clusTrackIDs_tag = meas1.tag_prop{curr_file_index}.trackIDs(meas1.tag_prop{curr_file_index}.N>=10);
    curr_clusTrackIDs_bin = meas1.bin_prop{curr_file_index}.trackIDs(meas1.bin_prop{curr_file_index}.N>=10);

    clus_to_use_tag = ~annotated_propdist.allanal.propdists.tag{i}.CLVgt2(:) & annotated_propdist.allanal.propdists.tag{i}.startAfter4Sec(:) & ~annotated_propdist.allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ~annotated_propdist.allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:);
    clus_to_use_bin = ~annotated_propdist.allanal.propdists.bin{i}.CLVgt2(:) & annotated_propdist.allanal.propdists.bin{i}.startAfter4Sec(:) & ~annotated_propdist.allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ~annotated_propdist.allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:);

    curr_clusTrackIDs_tag = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_tag(clus_to_use_tag),'uniformoutput',false));
    curr_clusTrackIDs_bin = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_bin(clus_to_use_bin),'uniformoutput',false));

    % all usable track IDs
    tag_track_IDs = load(sprintf('%s\\vbspt_source\\tag_%02d.mat',curr_filelocation,curr_file_index),'vbspt_trackIDs');
    bin_track_IDs = load(sprintf('%s\\vbspt_source\\bin_%02d.mat',curr_filelocation,curr_file_index),'vbspt_trackIDs');
    
    % all vbSPT annotations
    tag_track_diffstates = load(sprintf('%s\\vbspt_result/tag_%02d_HMManalysis_hidden2.mat',curr_filelocation,curr_file_index),'Wbest');
    bin_track_diffstates = load(sprintf('%s\\vbspt_result/bin_%02d_HMManalysis_hidden2.mat',curr_filelocation,curr_file_index),'Wbest');
    
    nonclus_tag_trackIDs = setdiff(tag_track_IDs.vbspt_trackIDs,curr_clusTrackIDs_tag);
    nonclus_bin_trackIDs = setdiff(bin_track_IDs.vbspt_trackIDs,curr_clusTrackIDs_bin);
    
    % Determine what the appropriate vbSPT track ID is for each [all track] ID.
    curr_vbsptID_clusTracks_tag = zeros(numel(curr_clusTrackIDs_tag),1);
    curr_vbsptID_clusTracks_bin = zeros(numel(curr_clusTrackIDs_bin),1);
    curr_vbsptID_nonclusTracks_tag = zeros(numel(nonclus_tag_trackIDs),1);
    curr_vbsptID_nonclusTracks_bin = zeros(numel(nonclus_bin_trackIDs),1);
    
    for j=1:numel(curr_clusTrackIDs_tag)
        curr_vbsptID_clusTracks_tag(j) = find(tag_track_IDs.vbspt_trackIDs==curr_clusTrackIDs_tag(j),1);
    end
    for j=1:numel(curr_clusTrackIDs_bin)
        curr_vbsptID_clusTracks_bin(j) = find(bin_track_IDs.vbspt_trackIDs==curr_clusTrackIDs_bin(j),1);
    end
    for j=1:numel(nonclus_tag_trackIDs)
        curr_vbsptID_nonclusTracks_tag(j) = find(tag_track_IDs.vbspt_trackIDs==nonclus_tag_trackIDs(j),1);
    end
    for j=1:numel(nonclus_bin_trackIDs)
        curr_vbsptID_nonclusTracks_bin(j) = find(bin_track_IDs.vbspt_trackIDs==nonclus_bin_trackIDs(j),1);
    end
    
    tag_clus_totSlowFrames(i) = sum(cellfun(@(x) sum(double(x==1)),tag_track_diffstates.Wbest.est2.sMaxP(curr_vbsptID_clusTracks_tag)));
    bin_clus_totSlowFrames(i) = sum(cellfun(@(x) sum(double(x==1)),bin_track_diffstates.Wbest.est2.sMaxP(curr_vbsptID_clusTracks_bin)));
    tag_nonclus_totSlowFrames(i) = sum(cellfun(@(x) sum(double(x==1)),tag_track_diffstates.Wbest.est2.sMaxP(curr_vbsptID_nonclusTracks_tag)));
    bin_nonclus_totSlowFrames(i) = sum(cellfun(@(x) sum(double(x==1)),bin_track_diffstates.Wbest.est2.sMaxP(curr_vbsptID_nonclusTracks_bin)));
    
    tag_clus_totFrames(i) = sum(cellfun(@(x) numel(x),tag_track_diffstates.Wbest.est2.sMaxP(curr_vbsptID_clusTracks_tag)));
    bin_clus_totFrames(i) = sum(cellfun(@(x) numel(x),bin_track_diffstates.Wbest.est2.sMaxP(curr_vbsptID_clusTracks_bin)));
    tag_nonclus_totFrames(i) = sum(cellfun(@(x) numel(x),tag_track_diffstates.Wbest.est2.sMaxP(curr_vbsptID_nonclusTracks_tag)));
    bin_nonclus_totFrames(i) = sum(cellfun(@(x) numel(x),bin_track_diffstates.Wbest.est2.sMaxP(curr_vbsptID_nonclusTracks_bin)));
end
% Get the total # tracks in clusters after filtering -- use the annotated_propdists data

save('track_slowFrac_time','tag_clus_totSlowFrames','bin_clus_totSlowFrames','tag_nonclus_totSlowFrames','bin_nonclus_totSlowFrames',...
                           'tag_clus_totFrames','bin_clus_totFrames','tag_nonclus_totFrames','bin_nonclus_totFrames');
end
load('track_slowFrac_time');
datamatrix = [tag_clus_totSlowFrames(:)./tag_clus_totFrames(:),...
              tag_nonclus_totSlowFrames(:)./tag_nonclus_totFrames(:),...
              bin_clus_totSlowFrames(:)./bin_clus_totFrames(:),...
              bin_nonclus_totSlowFrames(:)./bin_nonclus_totFrames(:)]*100;

figure('position',[345   700   240   418]);
hold on
colors=[51 204 255;255 255 255;255 102 102;255 255 255]/255;
for i=1:4
bar(i,mean(datamatrix(:,i)),'facecolor',colors(i,:));
errorbar(i,mean(datamatrix(:,i)),std(datamatrix(:,i)),'k');
end
set(gca,'xtick',1:4,'fontsize',14);
[~,p12]=ttest2(datamatrix(:,1),datamatrix(:,2));
[~,p34]=ttest2(datamatrix(:,3),datamatrix(:,4));
[~,p13]=ttest2(datamatrix(:,1),datamatrix(:,3));
[~,p24]=ttest2(datamatrix(:,2),datamatrix(:,4));

sigstar_v3({[1 2],[3 4],[1 3],[2 4]},...
        4*[p12 p34 p13 p24]);

set(gca,'xticklabel',{'tagSrc-CA, clus','tagSrc-CA, nonclus.',...
                      'Binder, clus','Binder, nonclus'});
set(gca,'xticklabelrotation',90)
ylabel('slow state (%)')
savefig(['percent_slowstate']);
print(gcf,'-dtiff',['percent_slowstate' '.tif'],'-r300');
print(gcf,'-dsvg','-painters',['percent_slowstate' '.svg'],'-r300');
          
% [~,~,stats]=anova1(datamatrix,[],'off');
% TKmat=multcompare(stats,'display','off'); % tukey-kramer correction for multiple comparisons
% [sigpair,anova_p]=format_stats_for_sigstar(TKmat);          
% 
% figure('position',[345   894   497   224]);
% hold on;
% plotSpread(datamatrix,'distributionColors',[51 204 255;0 0 0;255 102 102;0 0 0]/255);
% ylabel('time in slow state (%)');
% set(gca,'fontsize',14,'xticklabel',{'tagSrc, clus','tagSrc, nonclus','Binder, clus','Binder, nonclus'},...
%         'xticklabelrotation',15,'ylim',[0 100]);
% 
% % Render significance comparison bars
% sigstar_v3(sigpair,anova_p);

% figure('position',[705   912   303   216]);
% subplot(1,2,1); hold on;
% histogram(cell2mat(tag_cluslifetimes)*20.00,'binedges',(0:.0200:2)*1000,'facecolor',[51 102 51]/255,'edgecolor','none','normalization','pdf');
% %ylabel('probability density')
% histogram(cell2mat(tag_noncluslifetimes)*20.00,'binedges',(0:.0200:2)*1000,'edgecolor','k','displaystyle','stairs','normalization','pdf');
% xlabel('track lifetime (ms)')
% ylabel('fraction of tracks');
% legend('clustered tagSrc','non-clustered tagSrc');
% subplot(1,2,2); hold on;
% histogram(cell2mat(bin_cluslifetimes)*20.00,'binedges',(0:.0200:2)*1000,'facecolor',[255 51 51]/255,'edgecolor','none','normalization','pdf');
% %ylabel('probability density')
% histogram(cell2mat(bin_noncluslifetimes)*20.00,'binedges',(0:.0200:2)*1000,'edgecolor','k','displaystyle','stairs','normalization','pdf');
% xlabel('track lifetime (ms)')
% ylabel('fraction of tracks');
% legend('clustered Binder','non-clustered Binder');

end

function [sigpairs,pvals] = format_stats_for_sigstar(stats)
% Input is the result from a Tukey-Kramer-corrected one-way ANOVA.
sigIDX = find(stats(:,end)<0.05);
sigpairs = {};
pvals = [];
for i=1:numel(sigIDX)
   sigpairs = [sigpairs, stats(sigIDX(i),1:2)];
   pvals = [pvals, stats(sigIDX(i),end)];
end
end
