function count_codiff_in_clusters()
load('checked_codiff_result');
load('/Users/mikepab/Documents/MATLAB/MATLAB/bindertag/20190513_noad_stringent_reanalyze/annotated_propdists.mat')

% Use the checked codiffusion result to know which codiffusion datasources
% match up to which cluster analysis result files
%       - this may be best done by hand

% Use the track IDs in the codiffusion matrices, and cross-reference against
% the track IDs in cluster assignments
%       - for now, just count whether or not this happened (out of all tracks)
%       - three categories: 'invalid' track based on the roi i chose for the cell,
%                           'nonclustered', or 'clustered'.


% need to clarify with bei why there are 19 entries here, despite only 15 cells
% and which cells they correspond to "real" stuff since only 14 of the entries
% actually have codiffusion events recorded.
% --> i can manually check through myself first
codif_0109 = load('/Volumes/Groups/Hahn_Lab/Lab Members folders/Bei/Data-BT/Binder-tag/01092018 - MEF/2L4H - Dox100pg/coincidence_tracks_new.mat');
codif_0128 = load('/Volumes/Groups/Hahn_Lab/Lab Members folders/Bei/Data-BT/Binder-tag/01282018 - MEF/SLBH-100pg/coincidence_tracks_new.mat');
codif_0325 = load('/Volumes/Groups/Hahn_Lab/Lab Members folders/Bei/Data-BT/Binder-tag/03252018 - SLBH/DOX_100pg/coincidence_tracks_new.mat');
codif_0328 = load('/Volumes/Groups/Hahn_Lab/Lab Members folders/Bei/Data-BT/Binder-tag/03282018 - SLBH/100pg/coincidence_tracks_new.mat');
codif_0405 = load('/Volumes/Groups/Hahn_Lab/Lab Members folders/Bei/Data-BT/Binder-tag/04052018/SLBH-100pg/coincidence_tracks_new.mat');
codif_0617 = load('/Volumes/Groups/Hahn_Lab/Lab Members folders/Bei/Data-BT/Binder-tag/06172018 - SLBH/SLBH-100pg/coincidence_tracks_new.mat');

codif_0328_0405.coindidence_tracks_ch1 = [codif_0328.coindidence_tracks_ch1;codif_0405.coindidence_tracks_ch1];
codif_0328_0405.coindidence_tracks_ch2 = [codif_0328.coindidence_tracks_ch2;codif_0405.coindidence_tracks_ch2];

mm0109 = [1 2 3 4 5 6 7 8 10 11 12 13 15 16 17;...
          1 2 3 4 5 6 7 8 10 11 12 13 15 16 17]';
mm0128 = [1 2 3 7 8;...
          1 2 3 5 6]';
mm0325 = [1 2 3 4 5 6 7 8 9 10 11 12;...
          1 2 3 4 5 6 7 8 9 10 11 12]';
mm0328_0405 = [1 2 3 4 5 6 7 8 9 10 11 12 13 15 16 18 19 21 22 23 24 25 26 27 28 29 30;...
               1 2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28]';
mm0617 = [1 2 3 4 5 6 7 8 9 10;...
          1 2 3 4 5 6 7 8 9 10]';

anal03280405.propdists.tag = allanal.propdists.tag(1:28);
anal03280405.propdists.bin = allanal.propdists.bin(1:28);
anal0109.propdists.tag = allanal.propdists.tag(29:45);
anal0109.propdists.bin = allanal.propdists.bin(29:45);
anal0128.propdists.tag = allanal.propdists.tag(46:51);
anal0128.propdists.bin = allanal.propdists.bin(46:51);
anal0617.propdists.tag = allanal.propdists.tag(52:61);
anal0617.propdists.bin = allanal.propdists.bin(52:61);

      
[n_clus_0109,n_codif_0109,n_clus_codif_0109]=annotate_codiff_events(codif_0109,mm0109,...
                       '/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180628_noAd_0109_100pgmeasured_cluster_props.mat',...
                       anal0109);
[n_clus_0128,n_codif_0128,n_clus_codif_0128]=annotate_codiff_events(codif_0128,mm0128,...
                       '/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180628_noAd_0128_100pgmeasured_cluster_props.mat',...
                       anal0128);
[n_clus_0617,n_codif_0617,n_clus_codif_0617]=annotate_codiff_events(codif_0617,mm0617,...
                       '/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180628_noAd_0617_100pgmeasured_cluster_props.mat',...
                       anal0617);
[n_clus_0328_0405,n_codif_0328_0405,n_clus_codif_0328_0405]=annotate_codiff_events(codif_0328_0405,mm0328_0405,...
                       '/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180524_0328-0405_noAd_100pg_measured_cluster_props.mat',...
                       anal03280405);
nco.tag = [n_codif_0109.tag,n_codif_0128.tag,n_codif_0328_0405.tag,n_codif_0617.tag]';
nclco.tag = [n_clus_codif_0109.tag,n_clus_codif_0128.tag,n_clus_codif_0328_0405.tag,n_clus_codif_0617.tag]';
nco.bin = [n_codif_0109.bin,n_codif_0128.bin,n_codif_0328_0405.bin,n_codif_0617.bin]';
nclco.bin = [n_clus_codif_0109.bin,n_clus_codif_0128.bin,n_clus_codif_0328_0405.bin,n_clus_codif_0617.bin]';
nclco.both = [n_clus_codif_0109.both,n_clus_codif_0128.both,n_clus_codif_0328_0405.both,n_clus_codif_0617.both]';

sum(nclco.both)./sum(nco.tag)
save('counted_codiff_in_clusters');
end

function [n_clus,n_codif,n_clus_codif] = annotate_codiff_events(codif,matching_matrix,measpropfile,analprop)
load(measpropfile);
%load(analpropfile);

codiff_idxs = matching_matrix(:,1);
clusresult_idxs = matching_matrix(:,2);

n_clus_codif.tag = zeros(numel(codiff_idxs,1));
n_codif.tag = zeros(numel(codiff_idxs,1));
n_clus.tag = zeros(numel(codiff_idxs,1));

n_clus_codif.bin = zeros(numel(codiff_idxs,1));
n_codif.bin = zeros(numel(codiff_idxs,1));
n_clus.bin = zeros(numel(codiff_idxs,1));

n_clus_codif.both = zeros(numel(codiff_idxs,1));

%fprintf('Annotating codiffusion events for...\n\t%s\n\t%s\n',measpropfile,analpropfile);
fprintf('Annotating codiffusion events for...\n\t%s\n',measpropfile);
for i=1:numel(codiff_idxs) % the # of actual data points in the cluster analysis
    curr_cluster_cell_id = clusresult_idxs(i);
    curr_codif_cell_id = codiff_idxs(i);
    
    clusters_to_use_tag = ~analprop.propdists.tag{curr_cluster_cell_id}.CLVgt2(:) & ...
                        analprop.propdists.tag{curr_cluster_cell_id}.startAfter4Sec(:) & ...
                        ~analprop.propdists.tag{curr_cluster_cell_id}.allLoc.KDE.multiRegion50(:) & ...
                        ~analprop.propdists.tag{curr_cluster_cell_id}.allLoc.KDE.multiRegion95(:);  
                    
    clusters_to_use_bin = ~analprop.propdists.bin{curr_cluster_cell_id}.CLVgt2(:) & ...
                        analprop.propdists.bin{curr_cluster_cell_id}.startAfter4Sec(:) & ...
                        ~analprop.propdists.bin{curr_cluster_cell_id}.allLoc.KDE.multiRegion50(:) & ...
                        ~analprop.propdists.bin{curr_cluster_cell_id}.allLoc.KDE.multiRegion95(:);  
    
    tag_cluster_track_ids = cell2mat(cellfun(@(x) x(:),tag_prop{curr_cluster_cell_id}.trackIDs(clusters_to_use_tag),'uniformoutput',false));
    bin_cluster_track_ids = cell2mat(cellfun(@(x) x(:),bin_prop{curr_cluster_cell_id}.trackIDs(clusters_to_use_bin),'uniformoutput',false));
    
    tag_codif_track_ids = cellfun(@(x) x(1,end-1),codif.coindidence_tracks_ch2{curr_codif_cell_id});
    bin_codif_track_ids = cellfun(@(x) x(1,end-1),codif.coindidence_tracks_ch1{curr_codif_cell_id});
    
    
    
    n_codif.tag(i) = numel(tag_codif_track_ids);
    n_codif.bin(i) = numel(bin_codif_track_ids);
    
    n_clus.tag(i) = numel(tag_codif_track_ids);
    n_clus.bin(i) = numel(bin_codif_track_ids);
    
    n_clus_codif.tag(i) = numel(intersect(tag_codif_track_ids,tag_cluster_track_ids));
    n_clus_codif.bin(i) = numel(intersect(bin_codif_track_ids,bin_cluster_track_ids));
    
    which_tag_codif_in_clus = false(size(tag_codif_track_ids));
    which_bin_codif_in_clus = false(size(bin_codif_track_ids));
    
    temp_tag = intersect(tag_codif_track_ids,tag_cluster_track_ids);
    temp_bin = intersect(bin_codif_track_ids,bin_cluster_track_ids);
    for j=1:n_clus_codif.tag(i)
       which_tag_codif_in_clus(temp_tag(j) == tag_codif_track_ids) = true;
    end
    for j=1:n_clus_codif.bin(i)
       which_bin_codif_in_clus(temp_bin(j) == bin_codif_track_ids) = true;
    end
    
    n_clus_codif.both(i) = sum(which_tag_codif_in_clus|which_bin_codif_in_clus);
    
    
end

end
