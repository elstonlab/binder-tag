% COMPARE_VORONOICELL_SIZE_DISTRIBUTIONS_SIM50X
% Determines threshold densities for cluster segmentation.
%
% Input:
%   dateprefix, simfilesuffix, expfilesuffix: strings for file loading/naming, e.g.
%       '20180508_'
%       'sim_den_distr_50x'
%       'exp_den_distr'
%   threshSuffix: string for file + figure naming, e.g.
%       'voronoi_thresholds'
%   parentFigDir: string for parent figure directory, e.g.
%       'figures_20180508/'
%
% Saves a .mat file to the working directory:
%   [dateprefix threshSuffix '.mat']
%
% This sets thresholds for the segmentation step.
%
% Part of the cluster_segmentation.m pipeline.
function compare_voronoicell_size_distributions_sim50x(dateprefix,simfilesuffix,expfilesuffix,threshSuffix,parentFigDir)

tag_exp = load([dateprefix 'tag' expfilesuffix]);
bin_exp = load([dateprefix 'bin' expfilesuffix]);
tag_sim = load([dateprefix 'tag' simfilesuffix]);
bin_sim = load([dateprefix 'bin' simfilesuffix]);

figSubDir = 'vor_seg_clus_comp/';
parentFigDir = [parentFigDir figSubDir];
warning('off');
mkdir(parentFigDir);
warning('on');

% Get the localization areas in nm^2
tag_exp_areas = cellfun(@convert_density_to_areas,tag_exp.density_distrs,'uniformoutput',false);
bin_exp_areas = cellfun(@convert_density_to_areas,bin_exp.density_distrs,'uniformoutput',false);

ncells = numel(tag_exp_areas);
% the simulations are a nested cell array. There are 50 "cell array" elements
% per biological cell (ncells)
tag_sim_areas = cell(ncells,50);
bin_sim_areas = cell(ncells,50);
for i=1:ncells
    tag_sim_areas(i,:) = cellfun(@convert_density_to_areas,tag_sim.density_distrs{i},'uniformoutput',false);
    bin_sim_areas(i,:) = cellfun(@convert_density_to_areas,bin_sim.density_distrs{i},'uniformoutput',false);
end
tag_color = [0,0.4470,0.7410];
bin_color = [1,0,0];
sim_color = [0,0,0];

ncells = numel(tag_exp_areas);
tag_cutoffs = zeros(ncells,1);
bin_cutoffs = zeros(ncells,1);
binwidth = 100; % nm^2; corresponds to 10 nm pixel

% It's useful to check the single cell comparisons to make sure no problems
% occurred during threshold determination.
singleCell_comparison_name = [parentFigDir dateprefix threshSuffix '/singleCell_comparison'];
mkdir(singleCell_comparison_name);
for i=1:ncells
    tag_cutoffs(i) = get_cutoff_50x(tag_exp_areas{i},tag_sim_areas(i,:),binwidth);
    bin_cutoffs(i) = get_cutoff_50x(bin_exp_areas{i},bin_sim_areas(i,:),binwidth);
    plot_single_cell_smoothed_intersected_50x(tag_exp_areas{i},tag_sim_areas(i,:),tag_color,sim_color,'Tag-Exp','Tag-Sim',binwidth,[0 2000000],sprintf('%s_%s_%02d',singleCell_comparison_name,'_tag',i));
    plot_single_cell_smoothed_intersected_50x(bin_exp_areas{i},bin_sim_areas(i,:),bin_color,sim_color,'Bin-Exp','Bin-Sim',binwidth,[0 2000000],sprintf('%s_%s_%02d',singleCell_comparison_name,'_bin',i));
end

% Save these cutoffs
tag_cutoffs = convert_areas_to_density(tag_cutoffs);
bin_cutoffs = convert_areas_to_density(bin_cutoffs);

save([dateprefix threshSuffix],'tag_cutoffs','bin_cutoffs','binwidth');

% Check variability of the density threshold across biological cells.
% Outliers may suggest something went wrong (either with the thresholding,
% or with the data)
f=figure('position',[101   827   469   415]);
subplot(1,2,1); hold on
plotSpread(tag_cutoffs,'distributionMarker','o');
children=get(gca,'children');
children.MarkerFaceColor = children.Color;
boxplot(tag_cutoffs,'colors','k','symbol','k','whisker',0)
set(gca,'xticklabel','Tag','xtick',1);
ylabel('Density map threshold (localizations/nm^2)')
set(gca,'fontsize',14,'fontname','arial')
subplot(1,2,2); hold on
plotSpread(bin_cutoffs,'distributionMarker','o','distributionColor',[1 0 0]);
children=get(gca,'children');
children.MarkerFaceColor = children.Color;
boxplot(bin_cutoffs,'colors','k','symbol','k','whisker',0)
set(gca,'xticklabel','Binder','xtick',1);
ylabel('Density map threshold (localizations/nm^2)')
set(gca,'fontsize',14,'fontname','arial')

overall_comparison_name = [parentFigDir dateprefix threshSuffix '/overall_comparison'];

% Save the figure
savefig(f,[overall_comparison_name '.fig']);
print(f,'-dtiff',[overall_comparison_name '.tif'],'-r300');
close(f);
end

function cutoff = get_cutoff_50x(expdata,simdata50x,binwidth)
% Find the intesection of the probability distributions between the
% simulated mean and the experimental data. Apply smoothing prior to
% intersection checking; empirically chose 5pt Gaussian smooth.
% Check the single cell traces to verify that the cutoffs are consistent with
% the distributions' appearances.
[exp_distr,exp_edges]=histcounts(expdata,'binwidth',binwidth);
sim_distr=cell(50,1);
sim_edges=cell(50,1);

% Figure out the appropriate bin limits across the whole sim data
maxbin = min(cellfun(@max,simdata50x));
minbin = min(cellfun(@min,simdata50x));

for i=1:50
    [sim_distr{i},sim_edges]=histcounts(simdata50x{i},'binwidth',binwidth,'binlimits',[minbin maxbin]);
end
mean_sim_distr = mean(cell2mat(sim_distr));
exp_cents = get_bin_centers(exp_edges);
sim_cents = get_bin_centers(sim_edges);
exp_distr = smoothdata(exp_distr,'gaussian',5);
mean_sim_distr = smoothdata(mean_sim_distr,'gaussian',5);
[x_ints,~,~,~] = intersections(exp_cents,exp_distr,sim_cents,mean_sim_distr);
cutoff = x_ints(1);
end


function plot_single_cell_smoothed_intersected_50x(expdata,simdata50x,expcolor,simcolor,explabel,simlabel,binwidth,xplotwindow,savename)
f= figure('position',[356        1000         406         270]); hold on
[exp_distr,exp_edges]=histcounts(expdata,'binwidth',binwidth);
sim_distr=cell(50,1);
sim_edges=cell(50,1);

% Figure out the appropriate bin limits across the whole sim data
maxbin = min(cellfun(@max,simdata50x));
minbin = min(cellfun(@min,simdata50x));

for i=1:50
    [sim_distr{i},sim_edges]=histcounts(simdata50x{i},'binwidth',binwidth,'binlimits',[minbin maxbin]);
end
mean_sim_distr = mean(cell2mat(sim_distr));

exp_cents = get_bin_centers(exp_edges);
sim_cents = get_bin_centers(sim_edges);
exp_distr = smoothdata(exp_distr,'gaussian',5);
sim_distr = smoothdata(mean_sim_distr,'gaussian',5);
[x_ints,y_ints,~,~] = intersections(exp_cents,exp_distr,sim_cents,sim_distr);
hExp=plot(exp_cents,exp_distr,'linewidth',4,'color',expcolor);
hSim=plot(sim_cents,sim_distr,'linewidth',4,'color',simcolor);
plot(x_ints(1),y_ints(1),'mo','markerfacecolor','m','markersize',8)
legend([hExp,hSim],explabel,simlabel,'location','northeast');
xlabel('Area of Voronoi cell (nm^2)')
ylabel('Counts');
set(gca,'fontsize',16,'fontname','arial','xlim',xplotwindow);

savefig(f,[savename '.fig']);
print(f,'-dtiff',[savename '.tif'],'-r300');
close(f)
end

function centers = get_bin_centers(binedges)
centers = mean([binedges(1:end-1);binedges(2:end)]);
end

function areas = convert_density_to_areas(density)
% Convert from localizations/nm^2 to nm^2
areas=1./density;
end
function density = convert_areas_to_density(area)
% Convert from localizations/nm^2 to nm^2
density=1./area;
end
