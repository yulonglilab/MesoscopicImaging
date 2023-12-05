%% StimRspParcel_Ver1
% Running time:20220629
%   ======================================================================================
% Fei Deng,20220629,用于基于前期计算结果进一步计算不同状态下不同脑区的反应
% Fei Deng,20220629,用于光激活结果分析
% =======================================================================================
%   数据说明：
%   时间单位s
% Analysis by AccuSleep (1 = REM sleep, 2 = wakefulness, 3 = NREM sleep, 4 = undefined)
% =======================================================================================
close all;
% clearvars -except XX,clc;
clear,clc,

RspPathList = {'G:\WFFM1\20211001_TS60_Photoactivation\ResultsU_fitted405\5-HT3.0+jR_TS60_20211001_Photoactivation_FreqSeries_1\AllenMap\5-HT3.0+jR_TS60_20211001_Photoactivation_FreqSeries_1_G_BVcor_50Hz 10s_TrialAveRsp_CortexParcel.mat'};

TargetLabel = 'CortexParcel'; % 输入目标文件夹/文件特征字符串

[SaveFolder,SessionTitle,ext] = fileparts(RspPathList{1});
%%
RspMat = [];
% StatusLabels = [];
RspPath = RspPathList{1,1};
disp(['Importing: ',RspPath]);
%     load(RspPath);
load(RspPath,'CortexRegionLR');
RspMat = CortexRegionLR.Rsp;
IdxPeak = 99:103;
RspByState = mean(RspMat(:,IdxPeak),2,'omitnan');

Idx = strfind(SessionTitle,TargetLabel);
FileName = fullfile(SaveFolder,[SessionTitle(1:Idx-1),'StimParcel.mat']);
save(FileName,'RspByState','RspMat');
disp('Finished!');
