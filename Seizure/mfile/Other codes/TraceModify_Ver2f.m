%% TraceModify_Ver2f
% Running time:20220419
%   =======================================================================================
% Fei Deng,20211227
% 基于DeleteBreak修改
% Fei Deng,20220211,修改为自动识别多通道
% Fei Deng,20220313,Ratio改为Cor
% Fei Deng,20220313,增加条件判断modify类型
% Fei Deng,20220328,增加对MFI数据的整理
%   ========================================== =============================================
close all;
clear,clc,
ChNum = 2; % 激发光总数
PeriodNumRangeS = 1; % 根据实际需要手动修改实际分析图片在文件夹中的帧数范围，PeriodNumRangeS为起始帧的编号
BreakIdxList = [9760,9900;16180,16400;30000,30180];  % 中断帧的范围
BreakIdxList = floor(BreakIdxList/ChNum);
TargetLabel = 'RspOfSession.mat'; % 输入目标文件夹/文件特征字符串
MFIflag = 1; % 1-分析MFI，0-不分析MFI
ModifyFlag = 2; % 1-做binning,用于复制做origin图，2-(不做binning)用于Spike2汇总All trace文件

switch ModifyFlag
    case 1
        disp('Generate data for Origin');
        ModifyValue = NaN;
    case 2
        disp('Generate data for Spike2');
        %         BreakIdxList = [];  % 中断帧的范围
        ModifyValue = 1;
end
%%
MfileDir = pwd;
addpath(genpath(MfileDir));
ParentFolder = fileparts(MfileDir);
cd(ParentFolder);
[file,path] = uigetfile(['*',TargetLabel],'请选择 RspOfSession.mat文件'); %文件路径
ChTagIdx = strfind(file,TargetLabel)-2;

TargetFileName = file;
TargetFileName(ChTagIdx) = '*';
cd(path);
targetPack = dir(TargetFileName); %待整理文件列表获取
BreakNum = size(BreakIdxList,1);

for fi = 1:size(targetPack,1)
    FileName = targetPack(fi).name;
    ChTag = FileName(ChTagIdx);
    FilePathName = fullfile(path,FileName);
    load(FilePathName);
    for mi = 0:MFIflag
        switch mi
            case 0
                DataTag = 'Rsp';
                switch ModifyFlag
                    case 1
                        disp('Generate data for Origin');
                        Data = cat(2,RspOfSession.TimeInHrs,RspOfSession.Cor{1,1},RspOfSession.ChS{1,1},RspOfSession.ChNorm{1,1});
                        
                    case 2
                        disp('Generate data for Spike2');
                        Data = cat(2,RspOfSession.TimeInSec,RspOfSession.Cor{1,1},RspOfSession.ChS{1,1},RspOfSession.ChNorm{1,1});
                end
                
                IdxBad = find(abs(Data)==Inf);
                if IdxBad
                    disp('Inf replaced by 0 in: ');
                    disp(num2str(IdxBad));
                    Data(IdxBad) = 0;
                end
                if BreakNum
                    for i = 1:BreakNum
                        IdxStart = BreakIdxList(i,1)-PeriodNumRangeS;
                        IdxEnd = BreakIdxList(i,2)-PeriodNumRangeS+2;
                        ModifyValueF = mean(Data([IdxStart-59:IdxStart,IdxEnd:IdxEnd+59],2:end),1)*ModifyValue;
                        BreakTime(i) = mean(Data([IdxStart,IdxEnd],1));
                        IdxStart = IdxStart+1;
                        IdxEnd = IdxEnd-1;
                        Data(IdxStart:IdxEnd,2:end) = repmat(ModifyValueF,IdxEnd-IdxStart+1,1);
                    end
                end
                
                Xtime = RspOfSession.TimeInSec;
                figure,hold on;
                plot(Xtime,Data(:,3),'g'); plot(Xtime,Data(:,4),'b'); plot(Xtime,Data(:,2),'k');
                set(gca,'FontSize',10); % ; xlim([0 ImagingTimeM(end)])
                legend('raw','ref','cor');
                xlabel('Time (sec)','FontSize',12);
                ylabel('dF/F0','FontSize',12);
                title({'raw and fitted response';strrep(FileName(1:end-4),'_','\_')},'FontSize',12);
            case 1
                DataTag = 'MFI';
                if MFIflag
                    switch ModifyFlag
                        case 1
                            disp('Generate data for Origin');
                            Data = cat(2,RspOfSession.TimeInHrs,MFIOfSession.Cor{1,1},MFIOfSession.ChS{1,1},MFIOfSession.ChNorm{1,1});
                        case 2
                            disp('Generate data for Spike2');
                            Data = cat(2,RspOfSession.TimeInSec,MFIOfSession.Cor{1,1},MFIOfSession.ChS{1,1},MFIOfSession.ChNorm{1,1});
                    end
                    
                    IdxBad = find(abs(Data)==Inf);
                    if IdxBad
                        disp('Inf replaced by 0 in: ');
                        disp(num2str(IdxBad));
                        Data(IdxBad) = 0;
                    end
                    if BreakNum
                        for i = 1:BreakNum
                            IdxStart = BreakIdxList(i,1)-PeriodNumRangeS;
                            IdxEnd = BreakIdxList(i,2)-PeriodNumRangeS+2;
                            ModifyValueF = mean(Data([IdxStart-59:IdxStart,IdxEnd:IdxEnd+59],2:end),1)*ModifyValue;
                            BreakTime(i) = mean(Data([IdxStart,IdxEnd],1));
                            IdxStart = IdxStart+1;
                            IdxEnd = IdxEnd-1;
                            Data(IdxStart:IdxEnd,2:end) = repmat(ModifyValueF,IdxEnd-IdxStart+1,1);
                        end
                    end
                    
                    Xtime = RspOfSession.TimeInSec;
                    figure,hold on;
                    plot(Xtime,Data(:,3),'g'); plot(Xtime,Data(:,4),'b'); plot(Xtime,Data(:,2),'k');
                    set(gca,'FontSize',10); % ; xlim([0 ImagingTimeM(end)])
                    legend('raw','ref','cor');
                    xlabel('Time (sec)','FontSize',12);
                    ylabel('Fluorescence (au)','FontSize',12);
                    title({'raw and fitted MFI';strrep(FileName(1:end-4),'_','\_')},'FontSize',12);
                end
        end
        
        %         IdxBad = find(abs(Data)==Inf);
        %         if IdxBad
        %             disp('Inf replaced by 0 in: ');
        %             disp(num2str(IdxBad));
        %             Data(IdxBad) = 0;
        %         end
        %
        %         if BreakNum
        %             for i = 1:BreakNum
        %                 IdxStart = BreakIdxList(i,1)-PeriodNumRangeS+1;
        %                 IdxEnd = BreakIdxList(i,2)-PeriodNumRangeS+1;
        %                 Data(IdxStart:IdxEnd,:) = NaN;
        %                 BreakTime(i) = mean(Data([IdxStart-1,IdxEnd+1],1));
        %             end
        %         end
        BinStep = 5;
        BinNum = floor(size(Data,1)/BinStep);
        DataBin = zeros(BinNum,size(Data,2));
        for i = 1:BinNum
            IdxStart = (i-1)*BinStep+1;
            IdxEnd = i*BinStep;
            DataBin(i,:) = mean(Data(IdxStart:IdxEnd,:),1);
        end
        eval(['Data_',DataTag,'_',ChTag,' = Data;']);
        eval(['DataBin_',DataTag,'_',ChTag,' = DataBin;']);
    end
end
disp([FileName,': please copy data from Data and DataBin']);
cd(MfileDir);

