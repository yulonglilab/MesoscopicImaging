function IMseries = load_IMseries_from_datbin_file_Ver1(filename)
%   =======================================================================================
% Fei Deng,202202427, 用于读取以二进制存储的datbin文件中的多维矩阵
%   =======================================================================================
% input parameter is a file path of the IMseries*.databin file
[fp,fn,ext] = fileparts(filename);
dashPositions = find(fn=='_');
H = str2double(fn(dashPositions(1)+1:dashPositions(2)-1));
W = str2double(fn(dashPositions(2)+1:dashPositions(3)-1));
NF = str2double(fn(dashPositions(3)+1:dashPositions(4)-1));
NS = str2double(fn(dashPositions(4)+1:dashPositions(5)-1));
NT = str2double(fn(dashPositions(5)+1:dashPositions(6)-1));
%     dotPositions = find(fn=='.');
%     name = fn(dashPositions(7)+1:dotPositions-1);
f = fopen(filename, 'r');
IMseries=fread(f,Inf,'*uint16');
fclose(f);
IMseries = reshape(IMseries,H,W,NF,NS,NT);


