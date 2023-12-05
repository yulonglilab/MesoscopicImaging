function IMseries = load_IMseries_from_datbin_file(folder)
%   =======================================================================================
% Fei Deng,20220328, 用于读取以二进制存储的datbin文件中的多维矩阵
%   =======================================================================================
% input parameter is a string of directory contain the IMseries*.databin file
    fdir = dir([folder '\IMseries*.datbin']);
    if isempty(fdir)
        errordlg('Can not find IMseries*.datbin file!!!', 'File Error');
    else
        disp(fdir.name);
    end
    fn = fdir.name;
    dashPositions = find(fn=='_');
    H = str2double(fn(dashPositions(1)+1:dashPositions(2)-1)); 
    W = str2double(fn(dashPositions(2)+1:dashPositions(3)-1)); 
    NF = str2double(fn(dashPositions(3)+1:dashPositions(4)-1)); 
    NS = str2double(fn(dashPositions(4)+1:dashPositions(5)-1)); 
    NT = str2double(fn(dashPositions(5)+1:dashPositions(6)-1));
%     dotPositions = find(fn=='.');
%     name = fn(dashPositions(7)+1:dotPositions-1);
    fnfull = fullfile(fdir.folder,fdir.name);
    f = fopen(fnfull, 'r');
    IMseries=fread(f,Inf,'*uint16');
    fclose(f);
    IMseries = reshape(IMseries,H,W,NF,NS,NT);
    
    
