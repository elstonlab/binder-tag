f=dir('../runinput_files/*m');

for i = 1:numel(f)
    currfile = fullfile(f(i).folder, f(i).name);
    
    R=VB3_HMManalysis(currfile);
    %disp('Completed ',);
    disp(currfile)
end