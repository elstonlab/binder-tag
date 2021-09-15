function track_lifetime_in_clus()

filelocations = {'/Volumes/Groups/Elston Lab/MikePablo/20180716 - satadrive backup/20180524_bintag_0328-0405_mergecopy/inputdata/',...
                 '/Volumes/Groups/Elston Lab/MikePablo/20180716 - satadrive backup/20180628_more_noAd/to_LL/01092018_100pg/inputdata/',...
                 '/Volumes/Groups/Elston Lab/MikePablo/20180716 - satadrive backup/20180628_more_noAd/to_LL/01282018_100pg/inputdata/',...
                 '/Volumes/Groups/Elston Lab/MikePablo/20180716 - satadrive backup/20180628_more_noAd/to_LL/06172018_100pg/inputdata/'};

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

all_path_2.tag = {};
all_path_2.bin = {};
all_file_2.tag = {};
all_file_2.bin = {};
for i=1:1
    [curr_path,curr_file] = batchGetPath(filelocations{2},'RightCali_lvPALM','mat');
    all_path_2.tag = [all_path_2.tag;curr_path];
    all_file_2.tag = [all_file_2.tag;curr_file];
    [curr_path,curr_file] = batchGetPath(filelocations{2},'Left_lvPALM','mat');
    all_path_2.bin = [all_path_2.bin;curr_path];
    all_file_2.bin = [all_file_2.bin;curr_file];
end
assert(numel(all_path_2.tag)==numel(all_path_2.bin));

all_path_3.tag = {};
all_path_3.bin = {};
all_file_3.tag = {};
all_file_3.bin = {};
for i=1:1
    [curr_path,curr_file] = batchGetPath(filelocations{3},'RightCali_lvPALM','mat');
    all_path_3.tag = [all_path_3.tag;curr_path];
    all_file_3.tag = [all_file_3.tag;curr_file];
    [curr_path,curr_file] = batchGetPath(filelocations{3},'Left_lvPALM','mat');
    all_path_3.bin = [all_path_3.bin;curr_path];
    all_file_3.bin = [all_file_3.bin;curr_file];
end
assert(numel(all_path_3.tag)==numel(all_path_3.bin));

all_path_4.tag = {};
all_path_4.bin = {};
all_file_4.tag = {};
all_file_4.bin = {};
for i=1:1
    [curr_path,curr_file] = batchGetPath(filelocations{4},'RightCali_lvPALM','mat');
    all_path_4.tag = [all_path_4.tag;curr_path];
    all_file_4.tag = [all_file_4.tag;curr_file];
    [curr_path,curr_file] = batchGetPath(filelocations{4},'Left_lvPALM','mat');
    all_path_4.bin = [all_path_4.bin;curr_path];
    all_file_4.bin = [all_file_4.bin;curr_file];
end
assert(numel(all_path_4.tag)==numel(all_path_4.bin));

load('/Users/mikepab/Documents/MATLAB/MATLAB/bindertag/20190528_rectoclus/annotated_propdists_wRecToClus.mat')

N_tag = cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.tag,'uniformoutput',false);
N_bin = cellfun(@(x) x.N(~x.CLVgt2(:) & x.startAfter4Sec(:) & ~x.allLoc.KDE.multiRegion50(:) & ~x.allLoc.KDE.multiRegion95(:)),allanal.propdists.bin,'uniformoutput',false);

percent_tracks_in_clus_tag = zeros(61,1);
percent_tracks_in_clus_bin = zeros(61,1);

tag_cluslifetimes = cell(61,1);
bin_cluslifetimes = cell(61,1);
tag_noncluslifetimes = cell(61,1);
bin_noncluslifetimes = cell(61,1);

meas1 = load('/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180524_0328-0405_noAd_100pg_measured_cluster_props.mat');
meas2 = load('/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180628_noAd_0109_100pgmeasured_cluster_props.mat');
meas3 = load('/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180628_noAd_0128_100pgmeasured_cluster_props.mat');
meas4 = load('/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180628_noAd_0617_100pgmeasured_cluster_props.mat');
        

% For each cell
for i=1:61
    % # tracks in clusters after stringent filtering.
    curr_N_tag = sum(N_tag{i});
    curr_N_bin = sum(N_bin{i});
    
    % Get the total usable tracks (read the vbSPT source file and check the # of vbspt_trackIDs.
    if i>=1 && i<=28
        curr_filelocation = filelocations{1};
        curr_file_index = i;
        
        curr_clusTrackIDs_tag = meas1.tag_prop{curr_file_index}.trackIDs(meas1.tag_prop{curr_file_index}.N>=10);
        curr_clusTrackIDs_bin = meas1.bin_prop{curr_file_index}.trackIDs(meas1.bin_prop{curr_file_index}.N>=10);
        
        clus_to_use_tag = ~allanal.propdists.tag{i}.CLVgt2(:) & allanal.propdists.tag{i}.startAfter4Sec(:) & ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:);
        clus_to_use_bin = ~allanal.propdists.bin{i}.CLVgt2(:) & allanal.propdists.bin{i}.startAfter4Sec(:) & ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:);
        
        curr_clusTrackIDs_tag = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_tag(clus_to_use_tag),'uniformoutput',false));
        curr_clusTrackIDs_bin = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_bin(clus_to_use_bin),'uniformoutput',false));
        
        tag_expdata = load([all_path_1.tag{curr_file_index},'/',all_file_1.tag{curr_file_index}]);
        bin_expdata = load([all_path_1.bin{curr_file_index},'/',all_file_1.bin{curr_file_index}]);
        
        tag_cluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),tag_expdata.smLinked(curr_clusTrackIDs_tag));
        bin_cluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),bin_expdata.smLinked(curr_clusTrackIDs_bin));
        
    elseif i>28 && i<=45
        curr_filelocation = filelocations{2};
        curr_file_index = i-28;
        
        curr_clusTrackIDs_tag = meas2.tag_prop{curr_file_index}.trackIDs(meas2.tag_prop{curr_file_index}.N>=10);
        curr_clusTrackIDs_bin = meas2.bin_prop{curr_file_index}.trackIDs(meas2.bin_prop{curr_file_index}.N>=10);
                
        clus_to_use_tag = ~allanal.propdists.tag{i}.CLVgt2(:) & allanal.propdists.tag{i}.startAfter4Sec(:) & ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:);
        clus_to_use_bin = ~allanal.propdists.bin{i}.CLVgt2(:) & allanal.propdists.bin{i}.startAfter4Sec(:) & ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:);
        
        curr_clusTrackIDs_tag = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_tag(clus_to_use_tag),'uniformoutput',false));
        curr_clusTrackIDs_bin = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_bin(clus_to_use_bin),'uniformoutput',false));
        
        tag_expdata = load([all_path_2.tag{curr_file_index},'/',all_file_2.tag{curr_file_index}]);
        bin_expdata = load([all_path_2.bin{curr_file_index},'/',all_file_2.bin{curr_file_index}]);
        
        tag_cluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),tag_expdata.smLinked(curr_clusTrackIDs_tag));
        bin_cluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),bin_expdata.smLinked(curr_clusTrackIDs_bin));
        
    elseif i>45 && i<=51
        curr_filelocation = filelocations{3};
        curr_file_index = i-45;
        
        curr_clusTrackIDs_tag = meas3.tag_prop{curr_file_index}.trackIDs(meas3.tag_prop{curr_file_index}.N>=10);
        curr_clusTrackIDs_bin = meas3.bin_prop{curr_file_index}.trackIDs(meas3.bin_prop{curr_file_index}.N>=10);
                
        clus_to_use_tag = ~allanal.propdists.tag{i}.CLVgt2(:) & allanal.propdists.tag{i}.startAfter4Sec(:) & ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:);
        clus_to_use_bin = ~allanal.propdists.bin{i}.CLVgt2(:) & allanal.propdists.bin{i}.startAfter4Sec(:) & ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:);
        
        curr_clusTrackIDs_tag = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_tag(clus_to_use_tag),'uniformoutput',false));
        curr_clusTrackIDs_bin = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_bin(clus_to_use_bin),'uniformoutput',false));
        
        tag_expdata = load([all_path_3.tag{curr_file_index},'/',all_file_3.tag{curr_file_index}]);
        bin_expdata = load([all_path_3.bin{curr_file_index},'/',all_file_3.bin{curr_file_index}]);
        
        tag_cluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),tag_expdata.smLinked(curr_clusTrackIDs_tag));
        bin_cluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),bin_expdata.smLinked(curr_clusTrackIDs_bin));
        
    elseif i>51 && i<=61
        curr_filelocation = filelocations{4};   
        curr_file_index = i-51;     
        
        curr_clusTrackIDs_tag = meas4.tag_prop{curr_file_index}.trackIDs(meas4.tag_prop{curr_file_index}.N>=10);
        curr_clusTrackIDs_bin = meas4.bin_prop{curr_file_index}.trackIDs(meas4.bin_prop{curr_file_index}.N>=10);
        
                
        clus_to_use_tag = ~allanal.propdists.tag{i}.CLVgt2(:) & allanal.propdists.tag{i}.startAfter4Sec(:) & ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion50(:) & ~allanal.propdists.tag{i}.allLoc.KDE.multiRegion95(:);
        clus_to_use_bin = ~allanal.propdists.bin{i}.CLVgt2(:) & allanal.propdists.bin{i}.startAfter4Sec(:) & ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion50(:) & ~allanal.propdists.bin{i}.allLoc.KDE.multiRegion95(:);
        
        curr_clusTrackIDs_tag = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_tag(clus_to_use_tag),'uniformoutput',false));
        curr_clusTrackIDs_bin = cell2mat(cellfun(@(x) x(:),curr_clusTrackIDs_bin(clus_to_use_bin),'uniformoutput',false));
        
        tag_expdata = load([all_path_4.tag{curr_file_index},'/',all_file_4.tag{curr_file_index}]);
        bin_expdata = load([all_path_4.bin{curr_file_index},'/',all_file_4.bin{curr_file_index}]);
        
        tag_cluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),tag_expdata.smLinked(curr_clusTrackIDs_tag));
        bin_cluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),bin_expdata.smLinked(curr_clusTrackIDs_bin));
    end
    % all usable track IDs
    tag_track_IDs = load(sprintf('%s/vbspt/vbspt_source/tag_%02d.mat',curr_filelocation,curr_file_index),'vbspt_trackIDs');
    bin_track_IDs = load(sprintf('%s/vbspt/vbspt_source/bin_%02d.mat',curr_filelocation,curr_file_index),'vbspt_trackIDs');
    
    N_usable_tag = numel(tag_track_IDs.vbspt_trackIDs);
    N_usable_bin = numel(bin_track_IDs.vbspt_trackIDs);
    
    nonclus_tag_trackIDs = setdiff(tag_track_IDs.vbspt_trackIDs,curr_clusTrackIDs_tag);
    nonclus_bin_trackIDs = setdiff(bin_track_IDs.vbspt_trackIDs,curr_clusTrackIDs_bin);
    
    tag_noncluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),tag_expdata.smLinked(nonclus_tag_trackIDs));
    bin_noncluslifetimes{i} = cellfun(@(x) x(end,end)-x(1,end),bin_expdata.smLinked(nonclus_bin_trackIDs));
    
    percent_tracks_in_clus_tag(i) = curr_N_tag/N_usable_tag*100;
    percent_tracks_in_clus_bin(i) = curr_N_bin/N_usable_bin*100;    
end
% Get the total # tracks in clusters after filtering -- use the annotated_propdists data

save('track_lifetime_analysis','tag_cluslifetimes','bin_cluslifetimes','tag_noncluslifetimes','bin_noncluslifetimes');

figure('position',[345   894   497   224]);
subplot(1,2,1); hold on;
histogram(cell2mat(tag_cluslifetimes)*.02,'binedges',(0:.0200:2),'facecolor',[51 102 51]/255,'edgecolor','none','normalization','probability');
%ylabel('probability density')
histogram(cell2mat(tag_noncluslifetimes)*.02,'binedges',(0:.0200:2),'edgecolor','k','displaystyle','stairs','normalization','probability');
xlabel('track lifetime (s)')
ylabel('fraction of tagSrc tracks');
legend('clustered','non-clustered');
set(gca,'fontsize',14);

subplot(1,2,2); hold on;
histogram(cell2mat(bin_cluslifetimes)*.02,'binedges',(0:.0200:2),'facecolor',[255 51 51]/255,'edgecolor','none','normalization','probability');
%ylabel('probability density')
histogram(cell2mat(bin_noncluslifetimes)*.02,'binedges',(0:.0200:2),'edgecolor','k','displaystyle','stairs','normalization','probability');
xlabel('track lifetime (s)')
ylabel('fraction of Binder tracks');
legend('clustered','non-clustered');
set(gca,'fontsize',14);

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

