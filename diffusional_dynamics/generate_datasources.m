% Converts datafiles to vbSPT-usable files.
function generate_datasources()
%% Get file structures (lvPALM files)
datafolders = {'/Volumes/SATAdrive/20180524_bintag_0328-0405_mergecopy/'};
all_path.tag = {};
all_path.bin = {};
all_file.tag = {};
all_file.bin = {};
for i=1:numel(datafolders)
    [curr_path,curr_file] = batchGetPath(datafolders{i},'RightCali_lvPALM','mat');
    all_path.tag = [all_path.tag;curr_path];
    all_file.tag = [all_file.tag;curr_file];
    [curr_path,curr_file] = batchGetPath(datafolders{i},'Left_lvPALM','mat');
    all_path.bin = [all_path.bin;curr_path];
    all_file.bin = [all_file.bin;curr_file];
end
assert(numel(all_path.tag)==numel(all_path.bin));

%ROIfile = '/Volumes/med/Hahn_Lab/Lab Members/Bei/Data Backup/Binder-tag/03282018 - SLBH/roi.roi';
ncells = numel(all_path.tag);

% Generate source files
for i=1:ncells
    create_vbSPT_source(fullfile(all_path.tag{i},all_file.tag{i}),...
                        fullfile(all_path.tag{i},'roi.roi'),...
                        'vbspt_source/',sprintf('tag_%02d',i));

    create_vbSPT_source(fullfile(all_path.bin{i},all_file.bin{i}),...
                        fullfile(all_path.bin{i},'roi.roi'),...
                        'vbspt_source/',sprintf('bin_%02d',i));
end

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

function create_vbSPT_source(datafilename,roifilename,savedirectory,savename)
% Rejects tracks outside of the cell ROI
% Interpolates missing frames in tracks
% Saves the track coordinates in a new file, along with the original source
% information, and an list of the track indices within the original source.
% (The # tracks here may be less than the # tracks in the source)
fprintf('Processing %s\n-->\t%s\n',datafilename,roifilename);
expdata = load(datafilename);
roidata = ReadImageJROI(roifilename);

% For rejecting tracks outside of the cell.
cell_ps = polyshape(roidata.mnCoordinates(:,1),roidata.mnCoordinates(:,2));
vbspt_tracks = cell(numel(expdata.smLinked),1);
vbspt_trackIDs = zeros(numel(expdata.smLinked),1);
for i=1:numel(expdata.smLinked)
    % Reject track if its centroid is out of the cell.
    centroid = mean(expdata.smLinked{i}(:,1:2));
    if ~isinterior(cell_ps,centroid(1),centroid(2))
        continue;
    end
    
    allfeatureNew = func_fill_missing_pos(expdata.smLinked{i}, 1);
    vbspt_tracks{i} = allfeatureNew(:,1:2)*0.1067; % converting px to um
    vbspt_trackIDs(i) = i;
end
vbspt_trackIDs = vbspt_trackIDs(cellfun(@(x) ~isempty(x),vbspt_tracks));
vbspt_tracks = vbspt_tracks(cellfun(@(x) ~isempty(x),vbspt_tracks));
mkdir(savedirectory)
save([savedirectory savename],'vbspt_trackIDs','vbspt_tracks','datafilename','roifilename');
end