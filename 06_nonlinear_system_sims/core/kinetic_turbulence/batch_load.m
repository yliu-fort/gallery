%%load batch
clear; clc; close all;

d = dir;
folder = [];
for i= 1 : numel(d)
    
    if(startsWith(d(i).name,"Sun") == 1)
        folder = [folder {d(i).name}];
    end
end

%% read params
paramlist = [];
for i= 1 : numel(folder)
    fid = fopen([folder{i} '/_config']);
    C = textscan(fid, '%s\t%f\t%f\t%f','Headerlines', 1);
    
    param.OUTPUT_FOLDER = folder{i};
    for j = 1:numel(C{1})
        try
            param.(C{1}{j}) = C{2}(j);
        end
    end
    
    paramlist = [paramlist;param];
    
    fclose(fid);
end

clear C param fid i j 