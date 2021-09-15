function cluster_assoc_codiffFrac= count_codiff_in_clusters_negctrl(simulated_negctrl_file)
%load('neg_control_codiff_pairs');
load(simulated_negctrl_file);
load('/Users/mikepab/Documents/MATLAB/MATLAB/bindertag/20190513_noad_stringent_reanalyze/annotated_propdists.mat');


mm0109 = [1:15;...
          1 2 3 4 5 6 7 8 10 11 12 13 15 16 17]';
mm0128 = [1:5;...
          1 2 3 5 6]';
% mm0325 = [1 2 3 4 5 6 7 8 9 10 11 12;...
%           1 2 3 4 5 6 7 8 9 10 11 12]';
mm0328_0405 = [1:27;...
               1 2 3 4 5 6 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28]';
mm0617 = [1:10;...
          1 2 3 4 5 6 7 8 9 10]';
anal03280405.propdists.tag = allanal.propdists.tag(1:28);
anal03280405.propdists.bin = allanal.propdists.bin(1:28);
anal0109.propdists.tag = allanal.propdists.tag(29:45);
anal0109.propdists.bin = allanal.propdists.bin(29:45);
anal0128.propdists.tag = allanal.propdists.tag(46:51);
anal0128.propdists.bin = allanal.propdists.bin(46:51);
anal0617.propdists.tag = allanal.propdists.tag(52:61);
anal0617.propdists.bin = allanal.propdists.bin(52:61);

      
[n_clus_0109,n_codif_0109,n_clus_codif_0109]=annotate_codiff_events(neg_control_pairs.d0109,mm0109,...
                       '/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180628_noAd_0109_100pgmeasured_cluster_props.mat',...
                       anal0109);
[n_clus_0128,n_codif_0128,n_clus_codif_0128]=annotate_codiff_events(neg_control_pairs.d0128,mm0128,...
                       '/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180628_noAd_0128_100pgmeasured_cluster_props.mat',...
                       anal0128);
[n_clus_0617,n_codif_0617,n_clus_codif_0617]=annotate_codiff_events(neg_control_pairs.d0617,mm0617,...
                       '/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180628_noAd_0617_100pgmeasured_cluster_props.mat',...
                       anal0617);
[n_clus_0328_0405,n_codif_0328_0405,n_clus_codif_0328_0405]=annotate_codiff_events(neg_control_pairs.d0328_0405,mm0328_0405,...
                       '/Volumes/Groups/Elston Lab/MikePablo/20180504-bindertag/20190410_rerun_pipeline_and_binomanalysis/result/20180524_0328-0405_noAd_100pg_measured_cluster_props.mat',...
                       anal03280405);

nco.tag = [n_codif_0109.tag,n_codif_0128.tag,n_codif_0328_0405.tag,n_codif_0617.tag]';
nclco.tag = [n_clus_codif_0109.tag,n_clus_codif_0128.tag,n_clus_codif_0328_0405.tag,n_clus_codif_0617.tag]';
nco.bin = [n_codif_0109.bin,n_codif_0128.bin,n_codif_0328_0405.bin,n_codif_0617.bin]';
nclco.bin = [n_clus_codif_0109.bin,n_clus_codif_0128.bin,n_clus_codif_0328_0405.bin,n_clus_codif_0617.bin]';
nclco.both = [n_clus_codif_0109.both,n_clus_codif_0128.both,n_clus_codif_0328_0405.both,n_clus_codif_0617.both]';

cluster_assoc_codiffFrac=sum(nclco.both)./sum(nco.tag);

save('count_simulated_codiff_in_clus_negctrl');
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
   
    tag_codif_track_ids = codif.tag{curr_codif_cell_id};
    bin_codif_track_ids = codif.bin{curr_codif_cell_id};
    
    
    
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
