%% SleepWakeStateAnalysis_Ver1
% Running time:20220524
%   ======================================================================================
% Fei Deng,20220407,用于划分SleepWakeState
% =======================================================================================
%   数据说明：
%   时间单位s
% =======================================================================================
close all;
% clearvars -except XX,clc;
clear,clc,
% EEGfreq = 1000; % EEGfreq is 1000 Hz
EpochTime = 2.5; % EpochTime(s) EpochTime = 2.5;

PeriodRange =[0,99999999]; % 截取block的时间范围（s）

AccuSleepPath = 'I:\DF Codes\MATLAB code from other lab\AccuSleep-master';
addpath(genpath(AccuSleepPath)); % add searching path for AccuSleep files
trainedNetworkPath = 'I:\DF Codes\MATLAB codes\MATLAB code from other lab\AccuSleep-master\trainedNetworks';
trainedNetworkName = 'trainedNetwork2,5secEpochs.mat';
trainedNetworkFullName = fullfile(trainedNetworkPath,trainedNetworkName);
%% 根据EEG EMG划分睡眠状态
MfileDir = pwd;
addpath(genpath(MfileDir));
ParentFolder = fileparts(MfileDir);
cd(ParentFolder);
StimRcdFolder = uigetdir('','请选择指定的StimRcd文件夹');
cd(StimRcdFolder);
filetype = 'mat';
% [file,path] = uigetfile(['*.',filetype]); %文件路径
file = dir(['*.',filetype]); %文件路径
FileName = file.name;
load(FileName);
SaveFolder = fullfile(StimRcdFolder,'SleepState');
mkdir(SaveFolder);
%% SleepWake state analysis
IdxL = floor(PeriodRange(1)/EEG_DCrmv.interval);
IdxL = max(IdxL,1);
IdxR = ceil(PeriodRange(2)/EEG_DCrmv.interval);
IdxR = min(IdxR,EEG_DCrmv.length);
SessionTimeRange = [IdxL IdxR]*EEG_DCrmv.interval;
disp(['Data range (s): ',num2str(SessionTimeRange)]);
if exist('EEG_DCrmv','var')
    disp('With EEG DC remove');
    EEGraw = EEG_DCrmv;
else
    disp('Without EEG DC remove');
    EEGraw = EEG;
end
EEG = EEGraw.values(IdxL:IdxR);
EMGraw = EMG;
EMG = EMG.values(IdxL:IdxR);
save([SaveFolder,'\EEG_Accuformat','.mat'],'EEG');
save([SaveFolder,'\EMG_Accuformat','.mat'],'EMG');
FileCheck = dir('SleepStateLabel_AccuSleep.mat');
if isempty(FileCheck)
    labels = [];
    save([SaveFolder,'\SleepStateLabel_AccuSleep.mat'],'labels');
else
    load('SleepStateLabel_AccuSleep.mat');
end
FileCheck = dir(trainedNetworkName);
if isempty(FileCheck)
    disp(['Copy trainedNetworkName: ',trainedNetworkName]);
    trainedNetworkFullNameN = fullfile(SaveFolder,trainedNetworkName);
    copyfile(trainedNetworkFullName,trainedNetworkFullNameN);
end

FileCheck = dir('calibration.mat');
if isempty(FileCheck)
    calibrationData = [];
    save([SaveFolder,'\calibration.mat'],'calibrationData');
    %         AccuSleep_GUI_DF1;
    AccuSleep_GUI;
    pause,
else
    load('calibration.mat');
end
%%
if ~exist('SleepStateLabel','file')
    uiwait(msgbox('Please generate SleepStateLabel_AccuSleep.mat','Warning','modal'));
end
FileName = fullfile(SaveFolder,'SleepStateLabel_AccuSleep.mat');
load(FileName);

FileName = fullfile(SaveFolder,'SleepStateLabel.mat');
save(FileName,'labels','EpochTime');

disp('All Finished!');