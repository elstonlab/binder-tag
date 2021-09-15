% load in the usable track ids from the vbspt datafiles
% load in the codiffusion lists to see how many co-diffusion events to resample
% set rng to default at the beginning for reproducibility, and then use datasample
% apply codiffusion calculations as in count_codiff_in_clusters

function generate_codiff_negcontrol()
%create_trackinfo();
dat=load('trackid_trunc_info_for_negcontrol');
[n_codif] = get_n_codif();

n_codif = structfun(@(x) x.tag,n_codif,'uniformoutput',false);
trackids.d0109 = dat.trackids{1};
trackids.d0128 = dat.trackids{2};
trackids.d0617 = dat.trackids{3};
trackids.d0328_0405 = dat.trackids{4};

neg_control_pairs.d0109.tag = cell(numel(n_codif.d0109),1);
neg_control_pairs.d0109.bin = cell(numel(n_codif.d0109),1);
neg_control_pairs.d0128.tag = cell(numel(n_codif.d0128),1);
neg_control_pairs.d0128.bin = cell(numel(n_codif.d0128),1);
neg_control_pairs.d0617.tag = cell(numel(n_codif.d0617),1);
neg_control_pairs.d0617.bin = cell(numel(n_codif.d0617),1);
neg_control_pairs.d0328_0405.tag = cell(numel(n_codif.d0328_0405),1);
neg_control_pairs.d0328_0405.bin = cell(numel(n_codif.d0328_0405),1);

for realization_idx = 1:10
    rng(realization_idx);
    fields = {'d0109','d0128','d0617','d0328_0405'};
    for i=1:4
       currfield = fields{i};

       for j=1:numel(n_codif.(currfield))
           neg_control_pairs.(currfield).tag{j} = datasample(trackids.(currfield).tag{j},n_codif.(currfield)(j));
           neg_control_pairs.(currfield).bin{j} = datasample(trackids.(currfield).tag{j},n_codif.(currfield)(j));
       end
    end
    
    save(sprintf('neg_control_codiff_pairs_realiz%02d',realization_idx),'neg_control_pairs');
end
end

function create_trackinfo()
filepaths = {'/Volumes/Groups/Elston Lab/MikePablo/20180716 - satadrive backup/20180628_more_noAd/to_LL/01092018_100pg/inputdata/vbspt/vbspt_source/',...
             '/Volumes/Groups/Elston Lab/MikePablo/20180716 - satadrive backup/20180628_more_noAd/to_LL/01282018_100pg/inputdata/vbspt/vbspt_source/',...
             '/Volumes/Groups/Elston Lab/MikePablo/20180716 - satadrive backup/20180628_more_noAd/to_LL/06172018_100pg/inputdata/vbspt/vbspt_source/',...
             '/Volumes/Groups/Elston Lab/MikePablo/20180716 - satadrive backup/20180524_bintag_0328-0405_mergecopy/inputdata/vbspt/vbspt_source/'};
nfiles_per_path = [17,6,10,28];

% we need to drop out 0109/9, 0109/14, 0128/4, 0328-0405/14

trackids = cell(4,1);
for i=1:4
    [trackids{i}] = get_usable_track_ids(filepaths{i},nfiles_per_path(i));
end
save('trackid_info_for_negcontrol');
end

function [n_codif] = get_n_codif()
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

n_codif.d0109=get_ncodif(codif_0109,mm0109);
n_codif.d0128=get_ncodif(codif_0128,mm0128);
n_codif.d0617=get_ncodif(codif_0617,mm0617);
n_codif.d0328_0405=get_ncodif(codif_0328_0405,mm0328_0405);


end

function n_codif = get_ncodif(codif,matching_matrix)

codiff_idxs = matching_matrix(:,1);

n_codif.tag = zeros(numel(codiff_idxs,1));

n_codif.bin = zeros(numel(codiff_idxs,1));

for i=1:numel(codiff_idxs) % the # of actual data points in the cluster analysis
    curr_codif_cell_id = codiff_idxs(i);
    
    tag_codif_track_ids = cellfun(@(x) x(1,end-1),codif.coindidence_tracks_ch2{curr_codif_cell_id});
    bin_codif_track_ids = cellfun(@(x) x(1,end-1),codif.coindidence_tracks_ch1{curr_codif_cell_id});
    
    n_codif.tag(i) = numel(tag_codif_track_ids);
    n_codif.bin(i) = numel(bin_codif_track_ids);    
end

end


function [trackids] = get_usable_track_ids(path,ncells)
trackids.tag = cell(ncells,1);
trackids.bin = cell(ncells,1);
fprintf('Working on %s\n',path);
for i=1:ncells
    tagdat = load(sprintf('%s/tag_%02d',path,i),'vbspt_trackIDs');
    bindat = load(sprintf('%s/bin_%02d',path,i),'vbspt_trackIDs');
    trackids.tag{i} = tagdat.vbspt_trackIDs;
    trackids.bin{i} = bindat.vbspt_trackIDs;
end
end