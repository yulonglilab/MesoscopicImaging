%% SummarizeSheet_Ver1e
% Running time:20220515
%   =======================================================================================
% Fei Deng,20211011,used for Summarize Sheet from a cell type data
% SummarizeSheet_Ver1e,Fei Deng,20220515,可以读取双通道，按照不同需求整理数据
%   =======================================================================================
clear,clc;
TargetLabel = 'IMseriesOfSingleTrial.mat'; % 输入目标文件夹/文件特征字符串
MfileDir = pwd;
addpath(genpath(MfileDir));
ParentFolder = fileparts(MfileDir);
cd(ParentFolder);
[file,path] = uigetfile(['*',TargetLabel],'请选择 IMseriesOfSingleTrial.mat文件'); %文件路径
ChTagIdx = strfind(file,TargetLabel)-2;
TargetFileName = file;
TargetFileName(ChTagIdx) = '*';
cd(path);
targetPack = dir(TargetFileName); %待整理文件列表获取

AnalysisFlag = 5;
% 1: 汇总同一个ROI在不同刺激条件下的反应表格
% 2:汇总不同ROI在同一个刺激条件下的反应表格
% 3:汇总同一个ROI在不同刺激条件下的peak反应表格
% 4:汇总不同ROI在同一个刺激条件下的peak反应表格PeakValue_Mask_G
% 5:汇总同一个ROI在不同刺激条件下的peak（指定范围）反应表格
% 6:汇总不同ROI在同一个刺激条件下的peak（指定范围）反应表格
DataTag = 'Mask'; % ROI名称
ROIsID = 1; % 第几个ROI，ROI序号
StimID = 1; % 第几种刺激，刺激序号
PeakRange = 99:103; % 指定Peak的范围
%%
switch AnalysisFlag
    case 1
        %% 汇总同一个ROI在不同刺激条件下的反应表格
        disp('汇总同一个ROI在不同刺激条件下的反应表格...');
        for fi = 1:size(targetPack,1)
            FileName = targetPack(fi).name;
            ChTag = FileName(ChTagIdx);
            FilePathName = fullfile(path,FileName);
            load(FilePathName); % 读取RspOfSingleTrial变量，行列分别为StimModeNum和ROIsNum-1
            DataA = RspOfSingleTrial.Cor_M;
            %     DataA = RspOfSingleTrial.ChS_M;
            %             DataA = RspOfSingleTrial.ChNorm_M;
            [RowNum ColNum] = size(DataA);
            Data = [];
            for i = 1:RowNum
                Data = cat(2,Data,DataA{i,ROIsID}); % Average，SEM
            end
            eval(['Data_',DataTag,'_',ChTag,' = Data;']);
        end
        
    case 2
        %% 汇总不同ROI在同一个刺激条件下的反应表格
        disp('汇总同一个ROI在不同刺激条件下的反应表格...');
        for fi = 1:size(targetPack,1)
            FileName = targetPack(fi).name;
            ChTag = FileName(ChTagIdx);
            FilePathName = fullfile(path,FileName);
            load(FilePathName); % 读取RspOfSingleTrial变量，行列分别为StimModeNum和ROIsNum-1
            DataA = RspOfSingleTrial.Cor_M;
            %     DataA = RspOfSingleTrial.ChS_M;
            %     DataA = RspOfSingleTrial.ChNorm_M;
            [RowNum ColNum] = size(DataA);
            Data = [];
            for i = 1:ColNum
                Data = cat(2,Data,DataA{StimID,i}); % Average，SEM
            end
            eval(['Data_',DataTag,'_',ChTag,' = Data;']);
        end
        
    case 3
        %% 汇总同一个ROI在不同刺激条件下的peak反应表格
        disp('汇总同一个ROI在不同刺激条件下的peak反应表格...');
        for fi = 1:size(targetPack,1)
            FileName = targetPack(fi).name;
            ChTag = FileName(ChTagIdx);
            FilePathName = fullfile(path,FileName);
            load(FilePathName); % 读取RspOfSingleTrial变量，行列分别为StimModeNum和ROIsNum-1
            DataA = RspOfSingleTrial.Cor;
            %     DataA = RspOfSingleTrial.ChS;
            %     DataA = RspOfSingleTrial.ChNorm;
            [RowNum ColNum] = size(DataA);
            for i = 1:RowNum
                DataTemp = DataA{i,ROIsID};
                for ti  = 1:size(DataTemp,2)
                    [PeakValue(i,ti),PeakIdx(i,ti)] = max(movmean(DataTemp(:,ti),10));
                end
            end
            eval(['PeakValue_',DataTag,'_',ChTag,' = PeakValue;']);
            eval(['PeakIdx_',DataTag,'_',ChTag,' = PeakValue;']);
        end
        
    case 4
        %% 汇总不同ROI在同一个刺激条件下的peak反应表格
        disp('汇总不同ROI在同一个刺激条件下的peak反应表格...');
        for fi = 1:size(targetPack,1)
            FileName = targetPack(fi).name;
            ChTag = FileName(ChTagIdx);
            FilePathName = fullfile(path,FileName);
            load(FilePathName); % 读取RspOfSingleTrial变量，行列分别为StimModeNum和ROIsNum-1
            DataA = RspOfSingleTrial.Cor;
            %     DataA = RspOfSingleTrial.ChS;
            %     DataA = RspOfSingleTrial.ChNorm;
            [RowNum ColNum] = size(DataA);
            for i = 1:ColNum
                DataTemp = DataA{StimID,i};
                for ti  = 1:size(DataTemp,2)
                    [PeakValue(i,ti),PeakIdx(i,ti)] = max(movmean(DataTemp(:,ti),10));
                end
            end
            eval(['PeakValue_',DataTag,'_',ChTag,' = PeakValue;']);
            eval(['PeakIdx_',DataTag,'_',ChTag,' = PeakValue;']);
        end
        
    case 5
        %% 汇总同一个ROI在不同刺激条件下的peak（指定范围）反应表格
        disp('汇总同一个ROI在不同刺激条件下的peak（指定范围）反应表格...');
        for fi = 1:size(targetPack,1)
            FileName = targetPack(fi).name;
            ChTag = FileName(ChTagIdx);
            FilePathName = fullfile(path,FileName);
            load(FilePathName); % 读取RspOfSingleTrial变量，行列分别为StimModeNum和ROIsNum-1
            DataA = RspOfSingleTrial.Cor;
            %     DataA = RspOfSingleTrial.ChS;
            %     DataA = RspOfSingleTrial.ChNorm;
            [RowNum ColNum] = size(DataA);
            for i = 1:RowNum
                DataTemp = DataA{i,ROIsID};
                for ti  = 1:size(DataTemp,2)
                    PeakValue(i,ti) = mean(DataTemp(PeakRange,ti),1);
                end
            end
            eval(['PeakValue_',DataTag,'_',ChTag,' = PeakValue;']);
        end
        
    case 6
        %% 汇总不同ROI在同一个刺激条件下的peak（指定范围）反应表格
        disp('汇总不同ROI在同一个刺激条件下的peak（指定范围）反应表格...');
        for fi = 1:size(targetPack,1)
            FileName = targetPack(fi).name;
            ChTag = FileName(ChTagIdx);
            FilePathName = fullfile(path,FileName);
            load(FilePathName); % 读取RspOfSingleTrial变量，行列分别为StimModeNum和ROIsNum-1
            DataA = RspOfSingleTrial.Cor;
            %     DataA = RspOfSingleTrial.ChS;
            %     DataA = RspOfSingleTrial.ChNorm;
            [RowNum ColNum] = size(DataA);
            for i = 1:ColNum
                DataTemp = DataA{StimID,i};
                for ti  = 1:size(DataTemp,2)
                    PeakValue(i,ti) = mean(DataTemp(PeakRange,ti),1);
                end
            end
            eval(['PeakValue_',DataTag,'_',ChTag,' = PeakValue;']);
        end
        
end
cd(MfileDir);
disp('Finished!');
%% 20220117,used for Summarize Sheet from a matrix
% clear,clc;
% DataSheetNum =2;
% for i = 1:DataSheetNum
%     eval(['DataRaw',num2str(i),'=',num2str(i),';']);
% disp(['Please copy and paste raw data to DataRaw',num2str(i)]);
% end
% [RowNum ColNum] = size(DataRaw1);
% DataSheet = zeros(RowNum,ColNum*DataSheetNum);
% for i = 1:DataSheetNum
%     eval(['DataSheet(:,i:DataSheetNum:end) = DataRaw',num2str(i),';']);
% end
% disp('Finished!');