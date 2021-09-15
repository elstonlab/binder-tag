% Recursively search directories for datafiles of interest.

% All files follow this pattern:
% J:\Hahn_Users\Bei\000BinderTag\20201113\MEF-SrcYF\Exp20ms*\processing\*
%
% For example:
% J:\Hahn_Users\Bei\000BinderTag\20201113\MEF-SrcYF\Exp20ms-C1_2\processing\long_roi_lvPALM.mat
% J:\Hahn_Users\Bei\000BinderTag\20201113\MEF-SrcYF\Exp20ms-C1_2\processing\short_roi_lvPALM.mat
% J:\Hahn_Users\Bei\000BinderTag\20201113\MEF-SrcYF\Exp20ms-C1_2\processing\cell.roi
%
% 'short' denotes Binder.
% 'long' denotes tagSrc

function copy_target_data()
    tagSrc_files = rdir('J:\Hahn_Users\Bei\000BinderTag\20201113\MEF-SrcYF\*\processing\long*lvPALM.mat');
    Binder_files = rdir('J:\Hahn_Users\Bei\000BinderTag\20201113\MEF-SrcYF\*\processing\short*lvPALM.mat');
    roi_files = rdir('J:\Hahn_Users\Bei\000BinderTag\20201113\MEF-SrcYF\*\processing\cell.roi');

    
%     for i=1:numel(tagSrc_files)
%         srcpath = tagSrc_files(i).name;
%         %Get the cell name
%         s = regexp(srcpath, '\', 'split');
%         cellname = s{7}; % This is very hacky, but is fine since we shouldn't need to run this function much..
%         
%         destpath = cellname;
%         mkdir(destpath);
%         destpath = fullfile(destpath, 'RightCali_lvPALM.mat');
%         [success, msg, msgid] = copyfile(srcpath, destpath);
%     end
%     
%     for i=1:numel(Binder_files)
%         srcpath = Binder_files(i).name;
%         %Get the cell name
%         s = regexp(srcpath, '\', 'split');
%         cellname = s{7}; % This is very hacky, but is fine since we shouldn't need to run this function much..
%         
%         destpath = cellname;
%         mkdir(destpath);
%         destpath = fullfile(destpath, 'Left_lvPALM.mat');
%         [success, msg, msgid] = copyfile(srcpath, destpath);
%     end
    for i=1:numel(roi_files)
        srcpath = roi_files(i).name;
        % Get the cell name
        s = regexp(srcpath, '\', 'split');
        cellname = s{7}; % This is very hacky, but is fine since we shouldn't need to run this function much..
        
        destpath = cellname;
        mkdir(destpath);
        destpath = fullfile(destpath, 'roi.roi');
        [success, msg, msgid] = copyfile(srcpath, destpath);
    end
end
