function [allPathName, allFileName]= batchGetPath(rootfolder, expression, ext)
if ismac || isunix
    allTiffFiles = rdir([rootfolder '/**/*.' ext]);
elseif ispc 
    allTiffFiles = rdir([rootfolder '/**/*.' ext]);
end
j = 0;
for i = 1 : length(allTiffFiles)
    if regexpi(allTiffFiles(i).name, expression)
        j = j+1;
        allPathName{j, 1} = allTiffFiles(i).folder;
        if ismac || isunix
            tempName = split(allTiffFiles(i).name, '/');
        elseif ispc
            tempName = split(allTiffFiles(i).name, '\');
        end
        allFileName{j, 1} = tempName{end};
    end
end
end
