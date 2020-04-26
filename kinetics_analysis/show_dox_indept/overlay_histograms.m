function overlay_histograms()
% This function renders all the process times.

h1=cell(4,1);
h2=cell(4,1);
h3=cell(4,1);
h4=cell(4,1);
edge=cell(4,1);
cats=cell(4,1);

% Render using a smoothed kernel density estimate, or using the raw bar heights.
render_as_ksdensity=false;

[h1{1},h2{1},h3{1},h4{1},edge{1},cats{1}] = setup_250pg_data(render_as_ksdensity);
[h1{2},h2{2},h3{2},h4{2},edge{2},cats{2}] = setup_100pg_data(render_as_ksdensity);
[h1{3},h2{3},h3{3},h4{3},edge{3},cats{3}] = setup_10pg_data(render_as_ksdensity);
[h1{4},h2{4},h3{4},h4{4},edge{4},cats{4}] = setup_0pg_data(render_as_ksdensity);

if render_as_ksdensity
    figure;
    subplot(2,2,1);
    event_i = [h1{1};h1{2};h1{3};h1{4}]';
    plot(edge{1},event_i);
    legend('1x','4x','6x','10x')
    ylabel('probability density')
    set(gca,'fontsize',14,'xlim',[0 1500]);
    subplot(2,2,2);
    event_ii = [h2{1};h2{2};h2{3};h2{4}]';
    plot(edge{1},event_ii);
    legend('1x','4x','6x','10x')
    ylabel('probability density')
    set(gca,'fontsize',14,'xlim',[0 1500]);
    subplot(2,2,3);
    event_iii = [h3{1};h3{2};h3{3};h3{4}]';
    plot(edge{1},event_iii);
    legend('1x','4x','6x','10x')
    ylabel('probability density')
    set(gca,'fontsize',14,'xlim',[0 1500]);
    subplot(2,2,4);
    event_iv = [h4{1};h4{2};h4{3};h4{4}]';
    plot(edge{1},event_iv);
    legend('1x','4x','6x','10x')
    ylabel('probability density')
    set(gca,'fontsize',14,'xlim',[0 1500]);
else
    plot_edges = mean([edge{1}(1:end-1);edge{1}(2:end)]);


    figure('position',[605   969   408   301]);
    subplot(2,2,1);
    event_i = [h1{1};h1{2};h1{3};h1{4}]';
    plot(plot_edges,event_i);
    legend('1x','4x','6x','10x','fontsize',10)
    ylabel('fraction')
    set(gca,'fontsize',14,'xlim',[0 1500]);
    subplot(2,2,2);
    event_ii = [h2{1};h2{2};h2{3};h2{4}]';
    plot(plot_edges,event_ii);
    %legend('1x','4x','6x','10x','fontsize',10)
    set(gca,'fontsize',14,'xlim',[0 1500]);
    subplot(2,2,3);
    event_iii = [h3{1};h3{2};h3{3};h3{4}]';
    plot(plot_edges,event_iii);
    %legend('1x','4x','6x','10x','fontsize',10)
    ylabel('fraction')
    set(gca,'fontsize',14,'xlim',[0 1500]);
    xlabel('time (ms)');
    subplot(2,2,4);
    event_iv = [h4{1};h4{2};h4{3};h4{4}]';
    plot(plot_edges,event_iv);
    %legend('1x','4x','6x','10x','fontsize',10)
    set(gca,'fontsize',14,'xlim',[0 1500],'ylim',[0 0.2]);
    xlabel('time (ms)');

    savefig(['process_timings']);
    print(gcf,'-dtiff',['process_timings' '.tif'],'-r300');
    print(gcf,'-dsvg','-painters',['process_timings' '.svg'],'-r300');
end

figure('position',[421   969   183   301]);
bar(cell2mat(cats)')
xlabel('codiffusion category');
ylabel('percentage (%)');
set(gca,'fontsize',14)
legend('1x','4x','6x','10x','fontsize',10)
savefig(['codiff_categories']);
print(gcf,'-dtiff',['codiff_categories' '.tif'],'-r300');
print(gcf,'-dsvg','-painters',['codiff_categories' '.svg'],'-r300');

end

function [h1,h2,h3,h4,edge,cats] = setup_250pg_data(doKSestimate)

load('checked_codiff_result_250pg.mat');
[n_cat1_0109,n_cat2_0109,n_cat3_0109,n_cat4_0109] = get_codif_cats(codif_0109);
[n_cat1_0128,n_cat2_0128,n_cat3_0128,n_cat4_0128] = get_codif_cats(codif_0128);
[n_cat1_0325,n_cat2_0325,n_cat3_0325,n_cat4_0325] = get_codif_cats(codif_0325);
[n_cat1_0617,n_cat2_0617,n_cat3_0617,n_cat4_0617] = get_codif_cats(codif_0617);

all_cat1 = n_cat1_0109+n_cat1_0128+n_cat1_0325+n_cat1_0617;
all_cat2 = n_cat2_0109+n_cat2_0128+n_cat2_0325+n_cat2_0617;
all_cat3 = n_cat3_0109+n_cat3_0128+n_cat3_0325+n_cat3_0617;
all_cat4 = n_cat4_0109+n_cat4_0128+n_cat4_0325+n_cat4_0617;

cats=[all_cat1,all_cat2,all_cat3,all_cat4]./(all_cat1+all_cat2+all_cat3+all_cat4)*100;

[t1_0109,t2_0109,t3_0109,t4_0109] = get_class_times_v2(codif_0109);
[t1_0128,t2_0128,t3_0128,t4_0128] = get_class_times_v2(codif_0128);
[t1_0325,t2_0325,t3_0325,t4_0325] = get_class_times_v2(codif_0325);
[t1_0617,t2_0617,t3_0617,t4_0617] = get_class_times_v2(codif_0617);

all_t1 = [t1_0109;t1_0128;t1_0325;t1_0617];
all_t2 = [t2_0109;t2_0128;t2_0325;t2_0617];
all_t3 = [t3_0109;t3_0128;t3_0325;t3_0617];
all_t4 = [t4_0109;t4_0128;t4_0325;t4_0617];

all_t1 = all_t1(~isnan(all_t1))*20;
all_t2 = all_t2(~isnan(all_t2))*20;
all_t3 = all_t3(~isnan(all_t3))*20;
all_t4_withzeros = all_t4(~(isnan(all_t4)))*20;
all_t4 = all_t4(~(isnan(all_t4) | (all_t4==0)))*20;

if doKSestimate
    [h1,edge]=ksdensity(all_t1,0:20:4500);
    [h2,~]=ksdensity(all_t2,0:20:4500);
    [h3,~]=ksdensity(all_t3,0:20:4500);
    [h4,~]=ksdensity(all_t4,0:20:4500);
else
    [h1,edge]=histcounts(all_t1,'binedges',0:20:4500);
    [h2,~]=histcounts(all_t2,'binedges',0:20:4500);
    [h3,~]=histcounts(all_t3,'binedges',0:20:4500);
    [h4,~]=histcounts(all_t4,'binedges',0:20:4500);
    h1=h1./sum(h1);
    h2=h2./sum(h2);
    h3=h3./sum(h3);
    h4=h4./sum(h4);
end

end

function [h1,h2,h3,h4,edge,cats] = setup_100pg_data(doKSestimate)

load('checked_codiff_result_100pg.mat');
[n_cat1_0109,n_cat2_0109,n_cat3_0109,n_cat4_0109] = get_codif_cats(codif_0109);
[n_cat1_0128,n_cat2_0128,n_cat3_0128,n_cat4_0128] = get_codif_cats(codif_0128);
[n_cat1_0325,n_cat2_0325,n_cat3_0325,n_cat4_0325] = get_codif_cats(codif_0325);
[n_cat1_0328,n_cat2_0328,n_cat3_0328,n_cat4_0328] = get_codif_cats(codif_0328);
[n_cat1_0405,n_cat2_0405,n_cat3_0405,n_cat4_0405] = get_codif_cats(codif_0405);
[n_cat1_0617,n_cat2_0617,n_cat3_0617,n_cat4_0617] = get_codif_cats(codif_0617);

all_cat1 = n_cat1_0109+n_cat1_0128+n_cat1_0325+n_cat1_0328+n_cat1_0405+n_cat1_0617;
all_cat2 = n_cat2_0109+n_cat2_0128+n_cat2_0325+n_cat2_0328+n_cat2_0405+n_cat2_0617;
all_cat3 = n_cat3_0109+n_cat3_0128+n_cat3_0325+n_cat3_0328+n_cat3_0405+n_cat3_0617;
all_cat4 = n_cat4_0109+n_cat4_0128+n_cat4_0325+n_cat4_0328+n_cat4_0405+n_cat4_0617;

cats=[all_cat1,all_cat2,all_cat3,all_cat4]./(all_cat1+all_cat2+all_cat3+all_cat4)*100;

[t1_0109,t2_0109,t3_0109,t4_0109] = get_class_times_v2(codif_0109);
[t1_0128,t2_0128,t3_0128,t4_0128] = get_class_times_v2(codif_0128);
[t1_0325,t2_0325,t3_0325,t4_0325] = get_class_times_v2(codif_0325);
[t1_0328,t2_0328,t3_0328,t4_0328] = get_class_times_v2(codif_0328);
[t1_0405,t2_0405,t3_0405,t4_0405] = get_class_times_v2(codif_0405);
[t1_0617,t2_0617,t3_0617,t4_0617] = get_class_times_v2(codif_0617);

all_t1 = [t1_0109;t1_0128;t1_0325;t1_0328;t1_0405;t1_0617];
all_t2 = [t2_0109;t2_0128;t2_0325;t2_0328;t2_0405;t2_0617];
all_t3 = [t3_0109;t3_0128;t3_0325;t3_0328;t3_0405;t3_0617];
all_t4 = [t4_0109;t4_0128;t4_0325;t4_0328;t4_0405;t4_0617];

all_t1 = all_t1(~isnan(all_t1))*20;
all_t2 = all_t2(~isnan(all_t2))*20;
all_t3 = all_t3(~isnan(all_t3))*20;
all_t4_withzeros = all_t4(~(isnan(all_t4)))*20;
all_t4 = all_t4(~(isnan(all_t4) | (all_t4==0)))*20;

if doKSestimate
    [h1,edge]=ksdensity(all_t1,0:20:4500,'boundarycorrection','reflection');
    [h2,~]=ksdensity(all_t2,0:20:4500,'boundarycorrection','reflection');
    [h3,~]=ksdensity(all_t3,0:20:4500,'boundarycorrection','reflection');
    [h4,~]=ksdensity(all_t4,0:20:4500,'boundarycorrection','reflection');
else
    [h1,edge]=histcounts(all_t1,'binedges',0:20:4500);
    [h2,~]=histcounts(all_t2,'binedges',0:20:4500);
    [h3,~]=histcounts(all_t3,'binedges',0:20:4500);
    [h4,~]=histcounts(all_t4,'binedges',0:20:4500);
    h1=h1./sum(h1);
    h2=h2./sum(h2);
    h3=h3./sum(h3);
    h4=h4./sum(h4);
end

end


function [h1,h2,h3,h4,edge,cats] = setup_10pg_data(doKSestimate)

load('checked_codiff_result_10pg.mat');
[n_cat1_0109,n_cat2_0109,n_cat3_0109,n_cat4_0109] = get_codif_cats(codif_0109);
[n_cat1_0325,n_cat2_0325,n_cat3_0325,n_cat4_0325] = get_codif_cats(codif_0325);
[n_cat1_0617,n_cat2_0617,n_cat3_0617,n_cat4_0617] = get_codif_cats(codif_0617);

all_cat1 = n_cat1_0109+n_cat1_0325+n_cat1_0617;
all_cat2 = n_cat2_0109+n_cat2_0325+n_cat2_0617;
all_cat3 = n_cat3_0109+n_cat3_0325+n_cat3_0617;
all_cat4 = n_cat4_0109+n_cat4_0325+n_cat4_0617;

cats=[all_cat1,all_cat2,all_cat3,all_cat4]./(all_cat1+all_cat2+all_cat3+all_cat4)*100;

[t1_0109,t2_0109,t3_0109,t4_0109] = get_class_times_v2(codif_0109);
[t1_0325,t2_0325,t3_0325,t4_0325] = get_class_times_v2(codif_0325);
[t1_0617,t2_0617,t3_0617,t4_0617] = get_class_times_v2(codif_0617);

all_t1 = [t1_0109;t1_0325;t1_0617];
all_t2 = [t2_0109;t2_0325;t2_0617];
all_t3 = [t3_0109;t3_0325;t3_0617];
all_t4 = [t4_0109;t4_0325;t4_0617];

all_t1 = all_t1(~isnan(all_t1))*20;
all_t2 = all_t2(~isnan(all_t2))*20;
all_t3 = all_t3(~isnan(all_t3))*20;
all_t4_withzeros = all_t4(~(isnan(all_t4)))*20;
all_t4 = all_t4(~(isnan(all_t4) | (all_t4==0)))*20;

if doKSestimate
    [h1,edge]=ksdensity(all_t1,0:20:4500);
    [h2,~]=ksdensity(all_t2,0:20:4500);
    [h3,~]=ksdensity(all_t3,0:20:4500);
    [h4,~]=ksdensity(all_t4,0:20:4500);
else
    [h1,edge]=histcounts(all_t1,'binedges',0:20:4500);
    [h2,~]=histcounts(all_t2,'binedges',0:20:4500);
    [h3,~]=histcounts(all_t3,'binedges',0:20:4500);
    [h4,~]=histcounts(all_t4,'binedges',0:20:4500);
    h1=h1./sum(h1);
    h2=h2./sum(h2);
    h3=h3./sum(h3);
    h4=h4./sum(h4);
end

end

function [h1,h2,h3,h4,edge,cats] = setup_0pg_data(doKSestimate)

load('checked_codiff_result_0pg.mat');
[n_cat1_0109,n_cat2_0109,n_cat3_0109,n_cat4_0109] = get_codif_cats(codif_0109);
[n_cat1_0128,n_cat2_0128,n_cat3_0128,n_cat4_0128] = get_codif_cats(codif_0128);
[n_cat1_0325,n_cat2_0325,n_cat3_0325,n_cat4_0325] = get_codif_cats(codif_0325);
[n_cat1_0617,n_cat2_0617,n_cat3_0617,n_cat4_0617] = get_codif_cats(codif_0617);

all_cat1 = n_cat1_0109+n_cat1_0128+n_cat1_0325+n_cat1_0617;
all_cat2 = n_cat2_0109+n_cat2_0128+n_cat2_0325+n_cat2_0617;
all_cat3 = n_cat3_0109+n_cat3_0128+n_cat3_0325+n_cat3_0617;
all_cat4 = n_cat4_0109+n_cat4_0128+n_cat4_0325+n_cat4_0617;

cats=[all_cat1,all_cat2,all_cat3,all_cat4]./(all_cat1+all_cat2+all_cat3+all_cat4)*100;

[t1_0109,t2_0109,t3_0109,t4_0109] = get_class_times_v2(codif_0109);
[t1_0128,t2_0128,t3_0128,t4_0128] = get_class_times_v2(codif_0128);
[t1_0325,t2_0325,t3_0325,t4_0325] = get_class_times_v2(codif_0325);
[t1_0617,t2_0617,t3_0617,t4_0617] = get_class_times_v2(codif_0617);

all_t1 = [t1_0109;t1_0128;t1_0325;t1_0617];
all_t2 = [t2_0109;t2_0128;t2_0325;t2_0617];
all_t3 = [t3_0109;t3_0128;t3_0325;t3_0617];
all_t4 = [t4_0109;t4_0128;t4_0325;t4_0617];

all_t1 = all_t1(~isnan(all_t1))*20;
all_t2 = all_t2(~isnan(all_t2))*20;
all_t3 = all_t3(~isnan(all_t3))*20;
all_t4_withzeros = all_t4(~(isnan(all_t4)))*20;
all_t4 = all_t4(~(isnan(all_t4) | (all_t4==0)))*20;

if doKSestimate
    [h1,edge]=ksdensity(all_t1,0:20:4500);
    [h2,~]=ksdensity(all_t2,0:20:4500);
    [h3,~]=ksdensity(all_t3,0:20:4500);
    [h4,~]=ksdensity(all_t4,0:20:4500);
else
    [h1,edge]=histcounts(all_t1,'binedges',0:20:4500);
    [h2,~]=histcounts(all_t2,'binedges',0:20:4500);
    [h3,~]=histcounts(all_t3,'binedges',0:20:4500);
    [h4,~]=histcounts(all_t4,'binedges',0:20:4500);
    h1=h1./sum(h1);
    h2=h2./sum(h2);
    h3=h3./sum(h3);
    h4=h4./sum(h4);
end

end

function [t_class1,t_class2,t_class3,t_class4] = get_class_times_v2(codif)
codif_tag = codif.coindidence_tracks_ch2;
codif_bin = codif.coindidence_tracks_ch1;

ncells = numel(codif_tag);
n_max_events = sum(cellfun(@numel, codif_bin));
t_class1 = nan(n_max_events,1); % Duration of Binder/tagSrc co-diffusion if just Binder leaves while tagSrc still there.
t_class2 = nan(n_max_events,1); % Duration of Binder/tagSrc co-diffusion if both Binder and tagSrc leave simultaneously
t_class3 = nan(n_max_events,1); % Duration of tagSrc AFTER the class1 event happens.
t_class4 = nan(n_max_events,1); % Delay time (time between tagSrc appearing, and Binder associating).

% 20190725:
% The original class order for "get_class_times_v2" was:
% (1) Binder departs alone [inactivation]
% (2) Binder departs with tag [active Src PM dissoc]
% (3) tagSrc lifetime after Binder leaves [Src PM dissociation]
% (4) delay time between tagSrc and Binder arrival [activation]
%
% We want to re-order these to be consistent with other analyses, so now we
% have:
% (1) [activation]
% (2) [inactivation]
% (3) [Src PM dissociation]
% (4) [active Src PM dissociation]
%

counter=1;
for i=1:ncells
    if isempty(codif_bin{i})
       continue; % Nothing of interest here.
    end
    nevents = numel(codif_bin{i});

    for j=1:nevents
        tag_frames = codif_tag{i}{j}(:,end);
        bin_frames = codif_bin{i}{j}(:,end);
        start_of_codiff = max(codif_tag{i}{j}(1,end),codif_bin{i}{j}(1,end));

        % Binder and tagSrc left at the same time
        if abs(tag_frames(end) - bin_frames(end)) < 0.9 % error tolerance?
            t_class4(counter) = tag_frames(end) - start_of_codiff;
        else
            % Binder left before tagSrc
            if tag_frames(end) > bin_frames(end)
                t_class2(counter) = bin_frames(end) - start_of_codiff;
                t_class3(counter) = tag_frames(end) - bin_frames(end);
                % Binder arrived after tagSrc
                if tag_frames(1) < bin_frames(1)
                    t_class1(counter) = start_of_codiff - tag_frames(1);
                else % Binder arrive before/with tagSrc
                    t_class1(counter) = 0;
                end
            end
        end
        counter=counter+1;
    end
end
end

function [n_cat1,n_cat2,n_cat3,n_cat4] = get_codif_cats(codif)
codif_tag = codif.coindidence_tracks_ch2;
codif_bin = codif.coindidence_tracks_ch1;

ncells = numel(codif_tag);
n_cat1 = 0; % Duration of Binder/tagSrc co-diffusion if just Binder leaves while tagSrc still there.
n_cat2 = 0; % Duration of Binder/tagSrc co-diffusion if both Binder and tagSrc leave simultaneously
n_cat3 = 0; % Duration of tagSrc AFTER the class1 event happens.
n_cat4 = 0; % Delay time (time between tagSrc appearing, and Binder associating).

for i=1:ncells
    if isempty(codif_bin{i})
       continue; % Nothing of interest here.
    end

    nevents = numel(codif_bin{i});
    for j=1:nevents
        tag_frames = codif_tag{i}{j}(:,end);
        bin_frames = codif_bin{i}{j}(:,end);

        % Binder and tagSrc left at the same time
        if abs(tag_frames(end) - bin_frames(end)) < 0.9 % error tolerance?
            if tag_frames(1) < bin_frames(1) % Binder arrived after tagSrc
                n_cat2=n_cat2+1;
            else % Binder arrived with tagSrc
                n_cat4=n_cat4+1;
            end
        elseif tag_frames(end) > bin_frames(end) % Binder left before tagSrc
            if tag_frames(1) < bin_frames(1) % Binder arrived after tagSrc
                n_cat1=n_cat1+1;
            else % Binder arrived with tagSrc
                n_cat3=n_cat3+1;
            end
        end
    end
end
end
