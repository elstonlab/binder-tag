% Mike Pablo / 2018 Mar 9.
% ---------------------------------------------------------------------
% For each adhesion category, for tagSrc and Binder separately,
% measure the percentage of cells that enter adhesions directy via cytoplasm
% (first frame present in adhesion) or by lateral diffusion along the PM
% (first frame NOT present in adhesions);
% ---------------------------------------------------------------------
% v3 (Apr 24 2018) - switch order to [off na fc fa]. show ns bars
function percent_entry_memb_vs_cyto_noFAbound()
pre_calculated = false;
resultfilename = '20210328_percent_memb_vs_cyto_FAnoUpperBound';

if pre_calculated
    load(resultfilename);
else
    % Define/find the directories for the files of interest.
    datafolders = {'J:\Hahn_Users\Bei\Data-BT\0 new\For mike\20201125\',...
                   'J:\Hahn_Users\Bei\Data-BT\0 new\For mike\20201216\SLBH-Pxn\',...
                   'J:\Hahn_Users\Bei\Data-BT\0 new\For mike\20201219\SLBH-Pxn\'};
    all_path.tag = {};
    all_path.bin = {};
    all_file.tag = {};
    all_file.bin = {};
    
    % The first folder only has tagSrc data
    [curr_path,curr_file] = batchGetPath(datafolders{1},'corrDataStruct','mat');
    all_path.tag = [all_path.tag;curr_path];
    all_file.tag = [all_file.tag;curr_file];
    
    for i=2:numel(datafolders)
        [curr_path,curr_file] = batchGetPath(datafolders{i},'long_roi_scriptAnalysisResult_corrDataStruct','mat');
        all_path.tag = [all_path.tag;curr_path];
        all_file.tag = [all_file.tag;curr_file];
        [curr_path,curr_file] = batchGetPath(datafolders{i},'short_reg_roi_scriptAnalysisRe_corrDataStruct','mat');
        all_path.bin = [all_path.bin;curr_path];
        all_file.bin = [all_file.bin;curr_file];
    end
    

    ncells_tag = numel(all_path.tag);
    ncells_bin = numel(all_path.bin);
    
    % The percent_memb is 100-percent_cyto for each category.
    tag.NA.percent_cyto = zeros(ncells_tag,1);
    tag.FC.percent_cyto = zeros(ncells_tag,1);
    tag.FA.percent_cyto = zeros(ncells_tag,1);
    bin.NA.percent_cyto = zeros(ncells_bin,1);
    bin.FC.percent_cyto = zeros(ncells_bin,1);
    bin.FA.percent_cyto = zeros(ncells_bin,1);
    
    % Begin calculations
    for i=1:ncells_tag
       [tag.NA.percent_cyto(i),...
        tag.FC.percent_cyto(i),...
        tag.FA.percent_cyto(i)] = get_percent_cyto(fullfile(all_path.tag{i},all_file.tag{i}),...
                                                   fullfile(all_path.tag{i},'cell.roi'));
    end
    for i=1:ncells_bin
       [bin.NA.percent_cyto(i),...
        bin.FC.percent_cyto(i),...
        bin.FA.percent_cyto(i)] = get_percent_cyto(fullfile(all_path.bin{i},all_file.bin{i}),...
                                                   fullfile(all_path.bin{i},'cell.roi'));
    end

    % Save .mat files
    save(resultfilename);
end


% Clean up entries. Manual inspection showed that tagSrc cell 33
% had nan-valued recruitment and track density in nascent adhesion b/c 
% of no nascent adhesions present, so drop.
% For binder cell 16, all recruitment and track density entries were nan
% for adhesion categories, and Inf for non-adhesion categories. this was 
% because no adhesions were observed. We remove this row because
% all values are normalized to the average adhesion density.
cells_to_use_tag = [1:32, 34:38];
cells_to_use_bin = [1:15, 17:21];


render_stack_plots(tag, cells_to_use_tag,...
                   bin, cells_to_use_bin,...
                   'cyto_vs_memb_adhesion_entry');
                 
end

function render_stack_plots(tag, cells_to_use_tag, bin, cells_to_use_bin, figname)
%  Note that the input tag_density and bin_density are arranged in order of
% [NA FC FA nonad] so we permute the order first..
colors=[255,153,153;255,255,153;153,204,153]/255;
figure('position',[321 1057 267 235]);
subplot(1,2,1);
xvalues=[1 2 3 5 6 7];



for i=1:4
    bar(xvalues(i), mean(tag_density(:,i)), 'facecolor', colors(i,:)); hold on
    errorbar(xvalues(i), mean(tag_density(:,i)), nanstd(tag_density(:,i)), 'k');
end
ylabel(figname)
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','xtick',[],'ylim',[0 4]);
% Apply ANOVA to Tag and Binder groups separately, with Tukey-Kramer corr.
[~,~,tag_stats]=anova1(tag_density,[],'off');
TKmat_tag = multcompare(tag_stats,'display','off');
[sigpairs_tag,pvals_tag]=format_stats_for_sigstar(TKmat_tag);
sigstar_v3([sigpairs_tag],[pvals_tag]);

subplot(1,2,2);
xvalues=[1 2 3 4];
for i=1:4
    bar(xvalues(i), mean(bin_density(:,i)), 'facecolor', colors(i,:)); hold on
    errorbar(xvalues(i), mean(bin_density(:,i)), std(bin_density(:,i)), 'k');
end
ylabel(figname)
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','xtick',[],'ylim',[0 4]);
% Apply ANOVA to Tag and Binder groups separately, with Tukey-Kramer corr.
[~,~,tag_stats]=anova1(tag_density,[],'off');
TKmat_tag = multcompare(tag_stats,'display','off');
[sigpairs_tag,pvals_tag]=format_stats_for_sigstar(TKmat_tag);
sigstar_v3([sigpairs_tag],[pvals_tag]);



[~,~,bin_stats]=anova1(bin_density,[],'off');
TKmat_bin = multcompare(bin_stats,'display','off');
[sigpairs_bin,pvals_bin]=format_stats_for_sigstar(TKmat_bin);
%sigpairs_bin = cellfun(@(x) x+5,sigpairs_bin,'uniformoutput',false); % Align group index for plotting

% Render significance comparison bars
sigstar_v3([sigpairs_bin],[pvals_bin]);

% Update the legend at the very end to avoid automatic update weirdness
%lgd=legend('Off','NA','FC','FA','location','eastoutside');
%lgd.Box='off';
%mkdir(figdir)
box off;
savefig([figname]); % Save the .fig file
print(gcf,'-dtiff',[figname '.tiff'],'-r300'); % Save a 300 dpi tiff file
print(gcf,'-dpng',[figname '.png'],'-r300'); % Save a 300 dpi tiff file
print(gcf,'-painters','-dsvg',[figname '.svg'],'-r300'); % Save 300 dpi svg file
end

function render_flux_plots(tag_flux, bin_flux, figname)
colors=[255,153,153;255,255,153;153,204,153]/255;
figure('position',[321 1057 267 235]);
xvalues = [1 2 3 4 6 7 8 9];
subplot(1,2,1)
for i=1:3
    bar(xvalues(i),mean(tag_flux(:,i)),'facecolor',colors(i,:)); hold on
    errorbar(xvalues(i),mean(tag_flux(:,i)),std(tag_flux(:,i)),'k');
end
ylabel('memb. flux')
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','xtick',[],'ylim',[0 4]);
[~,~,tag_stats]=anova1(tag_flux,[],'off');
TKmat_tag = multcompare(tag_stats,'display','off');
[sigpairs_tag,pvals_tag]=format_stats_for_sigstar(TKmat_tag);

sigstar_v3([sigpairs_tag],[pvals_tag]);

subplot(1,2,2)
for i=1:3
    bar(xvalues(i),mean(bin_flux(:,i)),'facecolor',colors(i,:)); hold on
    errorbar(xvalues(i),mean(bin_flux(:,i)),std(bin_flux(:,i)),'k');
end
ylabel('memb. flux')
set(gca,'fontsize',14,'linewidth',1,'tickdir','out','xtick',[],'ylim',[0 4]);

% Apply ANOVA to Tag and Binder groups separately, with Tukey-Kramer corr.
[~,~,bin_stats]=anova1(bin_flux,[],'off');
TKmat_bin = multcompare(bin_stats,'display','off');
[sigpairs_bin,pvals_bin]=format_stats_for_sigstar(TKmat_bin);
% Render significance comparison bars
sigstar_v3([sigpairs_bin],[pvals_bin]);

box off;
savefig([figname]); % Save the .fig file
print(gcf,'-dtiff',[figname '.tiff'],'-r300'); % Save a 300 dpi tiff file
print(gcf,'-dpng',[figname '.png'],'-r300'); % Save a 300 dpi tiff file
print(gcf,'-painters','-dsvg',[figname '.svg'],'-r300'); % Save 300 dpi svg file
end

function [sigpairs,pvals] = format_stats_for_sigstar(stats)
% Input is the result from a Tukey-Kramer-corrected one-way ANOVA.
% sigIDX = find(stats(:,end)<0.05);
% sigpairs = {};
% pvals = [];
% for i=1:numel(sigIDX)
%    sigpairs = [sigpairs, stats(sigIDX(i),1:2)];
%    pvals = [pvals, stats(sigIDX(i),end)];
% end

sigpairs = cell(size(stats,1),1);
for i=1:size(stats,1)
    sigpairs{i} = stats(i,1:2);
end
pvals = stats(:,end);

end

function [allPathName, allFileName]= batchGetPath(rootfolder, expression, ext)
if ismac
    allTiffFiles = rdir([rootfolder '/**/*.' ext]);
elseif ispc 
    allTiffFiles = rdir([rootfolder '/**/*.' ext]);
end
j = 0;
for i = 1 : length(allTiffFiles)
    if regexpi(allTiffFiles(i).name, expression)
        j = j+1;
        allPathName{j, 1} = allTiffFiles(i).folder;
        if ismac
            tempName = split(allTiffFiles(i).name, '/');
        elseif ispc
            tempName = split(allTiffFiles(i).name, '\');
        end
        allFileName{j, 1} = tempName{end};
    end
end
end


function [PMentry,refdata] = get_perim_norm_PMrec(datafilename,roifilename)
% Returns 4x1 vectors of raw and normalized densities for the current cell.
% [NA FC FA nonad]; all normalized by the average adhesion density.
%
% refdata is a struct containing some useful extra information.
%   # total localizations per cell
%   # total tracks per cell
%   # total localizations per category per cell
%   # total tracks per category per cell

fprintf('Processing %s\n-->\t%s\n',datafilename,roifilename);
expdata = load(datafilename);
roidata = ReadImageJROI(roifilename);
cellarea = polyarea(roidata.mnCoordinates(:,1),roidata.mnCoordinates(:,2));

% For rejecting tracks outside of the cell.
cell_ps = polyshape(roidata.mnCoordinates(:,1),roidata.mnCoordinates(:,2));

nadframes = 1;

% NA FC FA nonad totalad
n_PMentry = zeros(5,1); 

area_ads = zeros(3,1);
perim_ads = zeros(3,1);
for i=1:numel(expdata.corrDataStruct)
    % Reject track if its centroid is out of the cell.
    centroid = mean(expdata.corrDataStruct{i}.position);
    if ~isinterior(cell_ps,centroid(1),centroid(2))
        continue;
    end
            
    % Check for tracks that diffuse into the adhesion from outside. 
    if (expdata.corrDataStruct{i}.flag_onAdhesion == 1) && (expdata.corrDataStruct{i}.d_min(1) > 0)
        adarea = expdata.focalAdhesion.Area(expdata.corrDataStruct{i}.adhesionID);
        %adperim = 
        if adarea >= 125 && adarea < 550
            n_PMentry(1) = n_PMentry(1)+1;
        end
        if adarea >= 550 && adarea < 2200
            n_PMentry(2) = n_PMentry(2)+1;
        end
        if adarea >= 2200 %&& adarea < 11000
            n_PMentry(3) = n_PMentry(3)+1;
        end
    elseif expdata.corrDataStruct{i}.flag_onAdhesion ~= 1
       n_PMentry(4) = n_PMentry(4)+1;
    end
end

n_PMentry(5) = sum(n_PMentry(1:3));

% Get adhesion areas
all_ad_area = 0;
all_ad_perim = 0;
for i=1:nadframes 
    for j=1:numel(expdata.focalAdhesion.Area)
        adarea = expdata.focalAdhesion.Area(j);
        adperim = perimeter(polyshape(expdata.focalAdhesion.Boundary{j}));
        
        if adarea >= 125 && adarea < 550
            area_ads(1) = area_ads(1)+adarea;
            perim_ads(1) = perim_ads(1)+adperim;
            all_ad_area = all_ad_area + adarea;
            all_ad_perim = all_ad_perim + adperim;
        end
        if adarea >= 550 && adarea < 2200
            area_ads(2) = area_ads(2)+adarea;
            perim_ads(2) = perim_ads(2)+adperim;
            all_ad_area = all_ad_area + adarea;
            all_ad_perim = all_ad_perim + adperim;
        end
        if adarea >= 2200 %&& adarea < 11000
            area_ads(3) = area_ads(3)+adarea;
            perim_ads(3) = perim_ads(3)+adperim;
            all_ad_area = all_ad_area + adarea;
            all_ad_perim = all_ad_perim + adperim;
        end
    end
end
% Get non-adhesion area
area_nonad = cellarea*25 - all_ad_area;
perim_nonad = nan; % This is not well-defined.

% Entry fluxes, arbitrary units.
PMentry = (n_PMentry ./ ([perim_ads; perim_nonad;sum(perim_ads)])) ./ (sum(n_PMentry(1:3))./sum(perim_ads));


% Set reference data
refdata.tot_tracks = numel(expdata.corrDataStruct);
refdata.tot_frames = sum(cellfun(@(x) size(x.position,1),expdata.corrDataStruct));

refdata.n_PMentry = n_PMentry;

refdata.area_ads = area_ads;
refdata.perim_ads = perim_ads;
refdata.nonad = area_nonad;
end


function [percent_cyto_NA, percent_cyto_FC, percent_cyto_FA] = get_percent_cyto(datafilename,roifilename)
% Returns 4x1 vectors of raw and normalized densities for the current cell.
% [NA FC FA nonad]; all normalized by the average adhesion density.
%
% refdata is a struct containing some useful extra information.
%   # total localizations per cell
%   # total tracks per cell
%   # total localizations per category per cell
%   # total tracks per category per cell

fprintf('Processing %s\n-->\t%s\n',datafilename,roifilename);
expdata = load(datafilename);
roidata = ReadImageJROI(roifilename);

% For rejecting tracks outside of the cell.
cell_ps = polyshape(roidata.mnCoordinates(:,1),roidata.mnCoordinates(:,2));

% NA FC FA
n_cyto_recruitment = zeros(3,1);
n_join_adhesion = zeros(3,1);
for i=1:numel(expdata.corrDataStruct)
    % Reject track if its centroid is out of the cell.
    centroid = mean(expdata.corrDataStruct{i}.position);
    if ~isinterior(cell_ps,centroid(1),centroid(2))
        continue;
    end
    
    % Check for tracks that start in ad.
    if expdata.corrDataStruct{i}.d_min(1) < 0
        adarea = expdata.focalAdhesion.Area(expdata.corrDataStruct{i}.adhesionID);
        if adarea >= 125 && adarea < 550
            n_cyto_recruitment(1) = n_cyto_recruitment(1)+1;
        end
        if adarea >= 550 && adarea < 2200
            n_cyto_recruitment(2) = n_cyto_recruitment(2)+1;
        end
        if adarea >= 2200 %&& adarea < 11000
            n_cyto_recruitment(3) = n_cyto_recruitment(3)+1;
        end
    end
        
    % Check for tracks that enter ad for >= 1 frame
    if expdata.corrDataStruct{i}.flag_onAdhesion == 1
        adarea = expdata.focalAdhesion.Area(expdata.corrDataStruct{i}.adhesionID);
        if adarea >= 125 && adarea < 550
            n_join_adhesion(1) = n_join_adhesion(1)+1;
        end
        if adarea >= 550 && adarea < 2200
            n_join_adhesion(2) = n_join_adhesion(2)+1;
        end
        if adarea >= 2200 %&& adarea < 11000
            n_join_adhesion(3) = n_join_adhesion(3)+1;
        end
    end
end
percent_cyto_NA = n_cyto_recruitment(1) / n_join_adhesion(1) * 100;
percent_cyto_FC = n_cyto_recruitment(2) / n_join_adhesion(2) * 100;
percent_cyto_FA = n_cyto_recruitment(3) / n_join_adhesion(3) * 100;
end
