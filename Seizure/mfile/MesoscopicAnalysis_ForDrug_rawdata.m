%% MesoscopicAnalysis_rawdata
% Running time:20230925
%   ======================================================================================
% StimRspAnalysis_Ver10,20230925，Fei Deng edited
% This Matlab script is used for analyzing  mesoscopic imaging data from a
% custom-made dual-color macroscope at Yulong Li lab
%   =======================================================================================
%%
close all,
clear,clc;

%% Parameter setting
% Please change parameters as noted according to the specific experiment
% (Parameters here are used for dual-color imaging of g5-HT3.0 and jRGECO1a during KA-induced seizure)

MotionCorFlag = 1; % whether do Motion correction, 1-Yes, 0-No
RegisterFlag = 1; % Whether register images from different cameras, 1-Yes, 0-No
BaselineCorFlag = 0; % Whether correct the baseline，1-use functions in Matlab toolbox，e.g. exp2; 2-use beads function in Matlab toolbox；3-use smooth function (very slow)
SpectraUnmixFlag = 1; % Whether do spectra unmixing, 1-Yes, 0-No
CorrectMethodFlag = 2; % Methods for correction based on 405 nm light excited images, 1-ratio; 2-fit then substract
WholeSessionRspFlag = 1; % Whether analyze response if the whole session, 1-Yes, 0-No
StimRspFlag = 0; % Whether analyze response according different stimuli
BasalFlagOfSession = 1;  % Methods for generating the baseline image for the whole session: 1-Specified frame；2-Based on the stimulus marker；3. The X% lower value for each pixel in the entire session; 4-The union set of baseline from the begin of session and each trial; 5-based on other labelled baseline, such as REM during sleep wake cycles
IdxBaseList = ones(20,2); % The range of the actual input and processed image frames for each session, initialized to the minimum value to avoid errors. The first column: starting frame, second column: ending frame
IdxBaseList(:,2) = 9999999999; % Initialized to the maximum value to avoid errors
IdxBaseList(1,:) = floor([1,35]); % Manually set the corresponding frame numbers of the baseline for the specified session (according to the image numbers of each channel)
forceSaveFlag = 1; % Whether to save IMseries file in binary format, 1-Yes, 0-No
SaveInSubFolderFlag = 0; % Whether to save the results of different sessions in separate folders, 1-Yes, 0-No
MFIrecordFlag = 1; % Whether to analyze the average fluorescence intensity of the entire session, 0 - not analyzed
SaveSingleChFlag = [1,1]; % Whether to save the response data of each camera (green, red) for each channel (uncorrected response data), 0 - not save, 1 - save
% mainPath = 'K:\Mouse in vivo\20210104_L21\SleepWake'; % Enter the range of folders to search
TargetLabel = 'Default'; % Enter the target folder/file feature string
MatchPattern = 1; % Whether scanDir function requires a full match of characters
TargetType = 1; % 1 - folder, 2 - file, 3 - folder and file
filetype = 'tif';% Image format
% #################### Imaging Information ####################
CameraNum = 2; % Number of cameras used
ChStatus = [1,1,1]; % 488nm, 561nm, 405nm for imaging
ChSim = [1,1,0]; % Channels for simultaneous imaging - 488nm, 561nm, 405nm
ChNorm = 3; % Normalized channel, 1-488nm, 2-561nm, 3-405nm
CameraTag = {'G','R','F'};
ChName = {'Ex488','Ex561','Ex405'};
ChNameTag = {'B','G','V'};  % Abbreviation of excitation light
ChColor = {'g','r','b'};
fixCh = 1; % Fixed channel number used for registration
ExposurePrecision = -2; % Adjusted based on the actual exposure time precision, e.g. 0.035s corresponds to precision -3, 0.04s corresponds to precision -2
ImSampleModeNum = 1; % Number of imaging sampling frequencies
Cgr = 0.119602778; % Ratio of green fluorescence signal crosstalk to red channel based on measured value, Cgr = 0.119602778 for g5-HT3.0
% Cgr = 0.098192171; % Ratio of green fluorescence signal crosstalk to red channel, based on measured and calibrated value, Cgr = 0.098192171 for eCB2.0
% Cgr = 0; % Ratio of green fluorescence signal crosstalk to red channel, based on measured and calibrated value, Cgr = 0.098192171 for eCB2.0; Theoretically, the green sensor will not cross-talked to red channel according to imaging setting
% Cgr = 0.073863753; % Ratio of green fluorescence signal crosstalk to red channel based on measured value, Cgr = 0.073863753 for EGFP
Crg = 0.001328723; % Ratio of red fluorescence signal crosstalk to green channel based on measured value, Cgr = 0.001328723 for jRGECO1a
% Crg = 2.69021E-07; % Ratio of red fluorescence signal crosstalk to green
% channel based on measured value,Cgr = 2.69021E-07 for r5-HT2.0
ChName = ChName(find(ChStatus==1));
ChColor = ChColor(find(ChStatus==1));
fixCh = sum(ChStatus(1:fixCh));
ChNumAll = length(ChName); % Total number or excitation light
ChSimNum = sum(ChSim);
ChNum = ChNumAll-(ChSimNum-1); % Number of excitation triggered per recording cycle
ChNorm0 = min([ChNum,ChNorm]); % 20230619 modify to ChNorm0 to distinguish it from ChNorm = 2 in the subsequent assignment and avoid errors caused by changes in the value
switch ChNum
    case 2
        ChTag0 = {'B','V';'G','V';};
    case 3
        ChTag0 = repmat(ChNameTag,CameraNum,1);
end

IMscale = 0.7; %Resize parameter of the image IMscale = 0.7;
PeriodNumRangeList = ones(10,2); % The range of actual input and processed image frames for each session, initialized to the minimum value to avoid errors
PeriodNumRangeList(:,2) = 9999999999; % Initialized to the maximum value to avoid errors
% PeriodNumRangeList(1,2) =200; % Manually modify PeriodNumRangeList according to actual needs
ImBinning = 4; % Binning value used during actual imaging ImBinning = 3;
pixelSzBin1 = 2211.358; % The scale during imaging while binning = 1X1 pixels / cm， calibrated in 20211031
pixelsize0 = pixelSzBin1/ImBinning;
pixelsize = pixelSzBin1/ImBinning*IMscale;
compression = 0;
% BgToAutoF = [2,2,2;2,2,2]; % Ratio of spontaneous fluorescence of the cortex and blood vessels in mice (no sensor expressed), each row represents a different camera channel, each column represents a different excitation light channel

% #################### Stimulation Information (can be left unchanged if no specific stimulation) ####################
TrialFastTime = 120;  % Duration of high-speed sampling for each trial (s)
StimT = 10;  % Maximum stimulation duration (s)
StimMaxT = StimT*1.2;  % Maximum stimulation duration (s)
StimTagList = {{'Base_20Hz 1s','Base_20Hz 10s','Base_50Hz 10s','GBR_20Hz 1s','GBR_20Hz 10s','GBR_50Hz 10s'}};% Stimulation name，one row per session
BaseDurList = ones(18,1)*10; %Time (s) used for baseline calculation for each trial of each session
IdxCorLen = 2; % To avoid the influence of photoactivation excitation light on the data used for baseline correction, the selected data points should shifted forward as much as possible. IdxCorLen is the number of data points shifted forward

% ResultTag = {'Miss','Hit','CR','FA'};

% #################### Fixed Parameters（Usually no change needed） ####################
BlockSz = 10000; % Process data in groups, set frames per group
MfileDir = pwd;
addpath(genpath(MfileDir));
ParentFolder = fileparts(MfileDir);
mfileNameR = mfilename; % Current .m file name
rawdataName = mfileNameR(strfind(mfileNameR,'rawdata'):end);
mainPath = fullfile(ParentFolder,rawdataName);
ImCorRatioFile = fullfile(MfileDir,'System Setting\CorRatio512pix_20211017.mat'); % Use this calibration file for data between 20211017 and 20220613
% ImCorRatioFile = fullfile(MfileDir,'CorRatio_20220613.mat'); % Use this calibration file for data after 20220613
load(ImCorRatioFile);
cd('I:\DF Codes\MATLAB codes\Motion_correction');
Mfiledir1 = pwd;
packDir1 = fileparts(Mfiledir1);
addpath(genpath(Mfiledir1)); %add path for searching functions
if CameraNum == 1
    RegisterFlag = 0;
end
switch CorrectMethodFlag
    case 1
        SaveFolder = fullfile(ParentFolder,['ResultsU',rawdataName(8:end)]);
    case 2
        SaveFolder = fullfile(ParentFolder,['ResultsU',rawdataName(8:end),'_fitted405']);
end
mkdir(SaveFolder);
%         cd('D:\Data\research\MATLAB codes\Motion_correction\NoRMCorre-master');

%% Scanning each session
tStart0 = tic;
cd(mainPath);
TargetList0 = scanDir(mainPath,TargetLabel,TargetType,MatchPattern);
TargetList0 = TargetList0';
fi = 0;
for i = 1:size(TargetList0,1)
    if isempty(strfind(TargetList0{i},'\Dark\'))
        fi = fi+1;
        TargetList(fi,1) = TargetList0(i);
    end
end
FileName = fullfile(SaveFolder,'TargetList.mat');
save(FileName);
TargetNum = size(TargetList,1);
SessionNum = TargetNum/CameraNum;

%% Loading shared image information
cd(TargetList{1});
targetPack = dir(['*.',filetype]);
IM1st = imread(targetPack(1).name);
[pixR,pixC] = size(IM1st); % Original image pixel size
% Read offset image
TargetPath = mainPath;
RefDarkA = uint16(zeros(pixR,pixC,CameraNum));
RefDarkA = F_RefDarkImport_Ver0(TargetPath,filetype,pixR,pixC,CameraNum,pixelSzBin1,ImBinning,compression);
%% Process different sessions
for fi = 1:SessionNum
    close all,
    TargetPath = TargetList{1+(fi-1)*CameraNum}; % All raw data is placed according to the session during recording
    %     TargetPath = TargetList{fi};
    [filepath,SessionTitle{fi,1},ext] = fileparts(fileparts(TargetPath));
    if ~isempty(ext) % Avoid errors caused by folder names containing decimal points
        SessionTitle{fi,1} = cat(2,SessionTitle{fi,1},ext);
    end
    TempTitle = SessionTitle{fi,1};
    Idx = find(TempTitle=='_',1,'last');
    SessionTitle{fi,1} = [TempTitle(1:Idx-3),TempTitle(Idx:end)];
    if SaveInSubFolderFlag
        SaveFolder1 = [SaveFolder,'\',SessionTitle{fi,1}];
        mkdir(SaveFolder1);
    else
        SaveFolder1 = SaveFolder;
    end
    SaveFolder2 = [SaveFolder1,'\Snapshot'];
    mkdir(SaveFolder2);
    % Read offset image
    TargetPathD = fileparts(fileparts(TargetPath));
    RefDark = F_RefDarkImport_Ver0(TargetPathD,filetype,pixR,pixC,CameraNum,pixelSzBin1,ImBinning,compression);
    if isempty(RefDark)
        if isempty(RefDarkA)
            disp(['No Dark images for ', SessionTitle{fi,1},'!']);
            RefDark = uint16(zeros(pixR,pixC,CameraNum));
        else
            RefDark = RefDarkA; % Use shared dark images
        end
    end
    %% Read Spike2 recording
    StimRcdPath = fullfile(filepath,'StimRcd'); % Path of StimRcd folder
    cd(StimRcdPath);
    Spk2RecordFile = dir('*.mat');
    Spk2RecordFile = Spk2RecordFile.name; % Spike2 recording task performance
    load(Spk2RecordFile);
    % ##########################################################################
    % Debug for extra point in spike2-recorded imaging data, delete the last extra point
    if  sfi == 999
        Imaging.level = Imaging.level(1:end-1);
        Imaging.times = Imaging.times(1:end-1);
    end
    % #####################################################################
    ImagingInfo = cat(2,Imaging.times,double(Imaging.level));
    % Imaging frequency calculation
    ImageStateLast = [ImagingInfo(2:end,1)-ImagingInfo(1:end-1,1)]; % Duration of each state, timing starts from the change of this state (0 or 1)
    ImageStateLastR = roundn(ImageStateLast,ExposurePrecision); %Adjust according to the actual exposure time precision
    %     ImageStateLastR = roundn(ImageStateLast,-4);
    ImageStateLastRFreq = tabulate(ImageStateLastR);
    ImageStateLastRFreq = sortrows(ImageStateLastRFreq,3,'descend');
    % For debugging when the number of frames recorded by spike2 is more than the actual imaging frames, find the extra data points as BugIdx1 to BugIdx21
    BugDurIdx = find(ImageStateLastRFreq(:,3)<0.01); % Time intervals with very low frequency should be abnormal
    while BugDurIdx
        if  BugDurIdx
            showboxh = errordlg(['Session ',num2str(fi),': error in Spike2 Imaging channel'],'Spike2 Error');
            waitfor(showboxh);
            BugFlag = 1;
            BugIdx = [ ];
            for bugi = 1:size(BugDurIdx,1)
                BugIdx = cat(1,BugIdx,find(ImageStateLastR==ImageStateLastRFreq(BugDurIdx(bugi),1)));
            end
            BugIdx = sort(BugIdx,'ascend');
            BugIdx = BugIdx+1; % Corresponding to the row number of ImagingInfo
        else
            BugFlag = 0;
        end
        if BugFlag
            figure,plot(ImageStateLastR(BugIdx(1)-1:max(BugIdx(end))+1));
            ImagingInfo0 = ImagingInfo;
            ImageStateLastR0 = ImageStateLastR;
            %         ImagingInfo = ImagingInfo0;
            %         figure,plot(ImagingInfo0(3000:3022,1),ImagingInfo0(3000:3022,2));
            ImagingInfo(BugIdx,:) = []; % For debugging when the number of frames recorded by spike2 is more than the actual imaging frames，delet extra data point
            %             ImagingInfo(BugIdx(1:end-1),:) = []; % For debugging when the number of frames recorded by spike2 is more than the actual imaging frames，delet extra data point
            %          ImagingInfo(BugIdx(1:end)+1,:) = []; % For debugging when spike2-recorded data is one frame more than the actual imaging frames，the extra data point is 66263，66264
            ImageStateLast = [ImagingInfo(2:end,1)-ImagingInfo(1:end-1,1)]; % Duration of each state, timing starts from the change of this state (0 or 1)
            ImageStateLastR = roundn(ImageStateLast,-2); %Adjust according to the actual exposure time precision
            ImageStateLast = [ImagingInfo(2:end,1)-ImagingInfo(1:end-1,1)]; % Duration of each state, timing starts from the change of this state (0 or 1)
            ImageStateLastR = roundn(ImageStateLast,-2); %Adjust according to the actual exposure time precision
            %     ImageStateLastR = roundn(ImageStateLast,-4);
            ImageStateLastRFreq = tabulate(ImageStateLastR);
            ImageStateLastRFreq = sortrows(ImageStateLastRFreq,3,'descend');
        end
        BugDurIdx = find(ImageStateLastRFreq(:,3)<0.01); % Time intervals with very low frequency should be abnormal
    end
    % Imaging Information
    ExTotalNum = size(ImagingInfo,1)/2; % Total number of excitation light pulses given
    ImagingTimeM = mean(reshape(ImagingInfo(:,1),ChNum*2,[])); % Calculate imaging time based on the number of channels
    BlockFrameNum = size(ImagingTimeM,2);
    %     ImSampleModeNum = size(ImageStateLastRFreq,1)-2;
    Period_T = zeros(1,ImSampleModeNum);
    if ChNum==1
        for i =  1:ImSampleModeNum
            Period_T(i) = ImageStateLastRFreq(1,1)+ImageStateLastRFreq(2+i,1); % Duration of each imaging cycle (including imaging and interval time), in seconds
        end
    elseif ChNum>1
        for i =  1:ImSampleModeNum
            Period_T(i) = ImageStateLastRFreq(1,1)*ChNum+ImageStateLastRFreq(2,1)*(ChNum-1)+ImageStateLastRFreq(2+i,1); % Duration of each imaging cycle (including imaging and interval time), in seconds
        end
    end
    Period_T = sort(Period_T,'ascend');
    Image_f = 1./Period_T; % Sampling frequency/Hz
    
    if StimRspFlag % Read stimulus information for analysis and organization
        StimInfo0 = cat(2,Stim.times,double(Stim.level));
        %         Stim.times(Stim.length-1:Stim.length) = []; % Debug for an extra high voltage record at the end
        %         Stim.level(Stim.length-1:Stim.length) = []; % Debug for an extra high voltage record at the end
        TagColor = {'black','red'}; % red-StimInfo, black-ImagingInfo
        % Calculate the number of stimuli and find the start point of each stimulus, generate StimInfo that only records the start and end time of each stimulus
        StimStateLast = [Stim.times(2:end)-Stim.times(1:end-1)]; % Duration of each state, timing starts from the change of this state (0 or 1)
        StimEndIdx = [find(StimStateLast>StimMaxT);size(Stim.times,1)]; % Index of the end time of each trial's stimulus in the original stimulus level record result
        TrialNum = size(StimEndIdx,1);
        StimStartIdx = [1;StimEndIdx(1:end-1)+1]; % Index of the start time of each trial's stimulus in the original stimulus level record result
        StimInfo = zeros(TrialNum*2,2);
        StimInfo(1:2:TrialNum*2-1,1) = Stim.times(StimStartIdx); % Start time of the stimulus in each trial
        StimInfo(1:2:TrialNum*2-1,2) = 1;
        %     StimInfo(2:2:TrialNum*2,1) = Stim.times(StimEndIdx); % End time of the stimulus in each trial
        MultiPulseIdx = find(StimEndIdx-StimStartIdx>1); % Find the trial numbers that contain multiple stimulus pulses
        if MultiPulseIdx % Determine if each trial has multiple stimulus pulses
            %              if size(StimInfo0,1)>size(StimInfo,1) % Determine if each trial has multiple stimulus pulses
            StimInfo(MultiPulseIdx'*2,1) = Stim.times(StimEndIdx(MultiPulseIdx))+(Stim.times(StimStartIdx(MultiPulseIdx)+2)-Stim.times(StimStartIdx(MultiPulseIdx)+1)); % At the end time of each trial's stimulus, add the time duration of non-high voltage pulses for correction
        else
            StimInfo(2:2:TrialNum*2,1) = Stim.times(StimEndIdx); % At the end time of each trial's stimulus, no correction needed when each trial has only one pulse
        end
        %         try
        %             StimInfo(2:2:TrialNum*2,1) = Stim.times(StimEndIdx)+(Stim.times(StimStartIdx+2)-Stim.times(StimStartIdx+1)); % At the end time of each trial's stimulus, add the time duration of non-high voltage pulses for correction
        %         catch exception
        %             StimInfo(2:2:TrialNum*2,1) = Stim.times(StimEndIdx); % At the end time of each trial's stimulus, no correction needed when each trial has only one pulse
        %         end
        StimInfo(2:2:TrialNum*2,2) = 0;
        
        %     if StimInfo(2,1)-StimInfo(1,1)>1.5*StimDur % If the time deviation of the first rectangular pulse is too large compared to the stimulus duration, it should be the pulse of the starting of arduino, discard it
        %         StimInfo = StimInfo(3:end,:);
        %     end
        
        %     TrialDur = 360;
        if TrialNum==1
            TrialDur = 360; % modify according to actual need？
        else
            TrialDur = (StimInfo(end-1,1)-StimInfo(1,1))/(TrialNum-1); %Average total time for each trial (s)
        end
        StimDur = StimInfo(2:2:end,1)-StimInfo(1:2:end-1,1); %Stimulus duration set for each trial, in seconds
        TrialDur = round(TrialDur); %The theoretical value is rounded to an integer
        StimDur = round(StimDur); %The theoretical value is rounded to an integer
        BaseDur = BaseDurList(fi,1); %Baseline time for each trial in each session set during analysis (s)
        RestDur = TrialDur-BaseDur; %Time from the start to the end of each stimulus (s)
        
        %     StimInfo = StimInfo(1:TrialNum*2,:); % Delete trials that were not recorded by arduino
        InfoData = {ImagingInfo,StimInfo};
        
        figure('Name','Stim & Imaging Info'),hold on;
        %     ylim([0 1.2]);
        TagColor = {'black','red'}; % red-StimInfo, black-ImagingInfo
        yFactor = [0.8,1];
        for i = 1:size(TagColor,2)
            InfoTemp = InfoData{i};
            xytmp = zeros(size(InfoTemp,1)*2,2);
            xytmp(2:2:end,:) = InfoTemp;
            xytmp(1,1) = InfoTemp(1,1);
            xytmp(3:2:end-1,1) =  xytmp(4:2:end,1);
            xytmp(3:2:end-1,2) =  xytmp(2:2:end-2,2);
            plot(xytmp(:,1),xytmp(:,2)*yFactor(i),'Color',TagColor{i},'LineWidth',1);
        end
        % Label the imaging moments
        yTemp = ones(1,BlockFrameNum)*0.85;
        text(ImagingTimeM,yTemp,num2cell(1:BlockFrameNum)); % Label the imaging period numbers
        StimTimeM = mean(reshape(StimInfo(:,1),2,[]));
        yTemp = ones(1,TrialNum)*0.95;
        text(StimTimeM,yTemp,num2cell(1:TrialNum),'Color','red'); % Label the imaging period numbers
        legend({'Imaging','Stim'});
        xlabel('Time (s)');
        saveas(gcf,[SaveFolder1,'\',SessionTitle{fi,1},'_Stim&Imaging Info.fig']);
        saveas(gcf,[SaveFolder1,'\',SessionTitle{fi,1},'_Stim&Imaging Info.jpg']);
        %% State labelling on images
        % The structure ImStamp records the labels and internal numbers of different events for each frame of the image, sorted by column, including Lick, Sound, Window, Actor, Result
        %     ImStamp.Flag = int8(zeros(BlockFrameNum,size(TaskData,2)));
        %     ImStamp.Note = int8(zeros(BlockFrameNum,size(TaskData,2)));
        ImStamp.Flag = zeros(BlockFrameNum,1)*nan;
        ImStamp.Note = ImStamp.Flag;
        % ####### Image labeling scheme 1: Attach each event label to the nearest frame of the event occurrence #######
        % Stimuli label
        ArdRcdFile = [Spk2RecordFile(1:end-4), '_ArdRcd.txt'];
        fullFileName = fullfile(StimRcdPath,ArdRcdFile);
        if exist(ArdRcdFile)
            SplitTag = 'StimMode: ';
            %             SplitTag = 'StimMode '; % Modify for abnormal recording file format
            %         SplitTag1 = 'StimInfo: ';
            %         SplitTag2 = ',';
            %                 StimMode = F_StimModeImport_Ver1(fullFileName,SplitTag1,SplitTag2);
            StimMode = F_StimModeImport_Ver0(fullFileName,SplitTag);
        else
            StimMode = ones(TrialNum,1);
        end
        StimModeNum = length(unique(StimMode));
        %         StimModeNum = max(StimMode);
        StimModeCount = zeros(1,StimModeNum);
        for i = 1:TrialNum
            CurImIdx = find(ImagingTimeM>StimInfo(2*i-1,1),1); % Find the image number that is greater than the start time of each stimulus
            if CurImIdx
                ImStamp.Flag(CurImIdx,1) = StimMode(i,1); % Event type number, i.e., stimulus mode number
                StimModeCount(StimMode(i,1)) = StimModeCount(StimMode(i,1))+1;
                ImStamp.Note(CurImIdx,1) = StimModeCount(StimMode(i,1)); %Event number, i.e., the number of the trial belonging to the event/stimulus
                ImStamp.Note(CurImIdx,2) = i; % Total trial number
            end
        end
        IdxStimAll = find(ImStamp.Note(:,2)>0);
        %     TrialFrames = round(mean(diff(IdxStimAll)));  % Average number of frames sampled per trial
        TrialFrames = min(diff(IdxStimAll));  % Minimum number of frames sampled per trial
        %         TrialFrames = max(diff(IdxStimAll));  % Maximum number of frames sampled per trial
        %         TrialFrames = round(mean(diff(IdxStimAll)));  % Average number of frames sampled per trial
        IdxStartAll = IdxStimAll-BaseDur*Image_f(1);
        IdxEndAll = IdxStartAll+TrialFrames-1;
        if length(Period_T)>1
            TrialFastFrames = TrialFastTime*Image_f(1); % Average number of frames of high-speed sampling per trial
            IdxFastSampleEndAll = IdxStartAll+TrialFastFrames-1;
            TrialTimeFrameRcd = cat(2,(1:TrialFastFrames)*Period_T(1),TrialFastFrames*Period_T(1)+(1:TrialFrames-TrialFastFrames)*Period_T(2)); % The corresponding recording time of each frame image for each trial
        else
            TrialTimeFrameRcd = (1:TrialFrames)*Period_T;
        end
        TrialTimeFrameRcd = TrialTimeFrameRcd-BaseDur;  % Set the stimulus moment as the zero moment, and the first frame after the stimulus as the first moment greater than zero
        TrialTimeFrameRcd = TrialTimeFrameRcd';
        save([SaveFolder1,'\',SessionTitle{fi,1},'_Stim&Imaging Info.mat'],'-v7.3');
    end
    
    %% Import imaging info
    ChAnalyse = zeros(CameraNum,ChNum); % 488nm, 561nm, 405nm for analyse
    for ci = 1:CameraNum % Process the data recorded in the same session separately by each camera ，CameraNum is the number of cameras
        switch ChSimNum
            case 1
                ChAnalyse(ci,[ci,ChNorm0]) = 1; % 488nm, 561nm, 405nm for analyse, for camera ci
            case 2
                ChAnalyse(ci,[1,ChNorm0]) = 1; % 488nm, 561nm, 405nm for analyse, for camera ci
        end
        ChAnlsNum(ci) = mean(sum(ChAnalyse(ci,:),2)); % channel number for analyse
        TargetPath = TargetList{ci+(fi-1)*CameraNum}; % The folder name of each camera recording is in the format of XXX_G_1, and all the raw data is placed according to the session
        %           TargetPath = TargetList{fi+(ci-1)*SessionNum};
        [filepath,SessionChTitle{fi,ci},ext] = fileparts(fileparts(TargetPath));
        if ~isempty(ext) % Avoid errors caused by folder names containing decimal points
            SessionChTitle{fi,ci} = cat(2,SessionChTitle{fi,ci},ext);
        end
        disp([num2str(fi),'#, Processing images in ',SessionChTitle{fi,ci}]),
        TempTitle = SessionChTitle{fi,ci};
        Idx = find(TempTitle=='_',1,'last');
        SessionChTitle{fi,ci} = [TempTitle(1:Idx-3),TempTitle(Idx:end),TempTitle(Idx-2:Idx-1)];
        ChTag(ci,:) = ChTag0(ci,find(ChAnalyse(ci,:)));
        cd(TargetPath);
        targetPack = dir(['*.',filetype]);
        targetsNameTemp = ({targetPack(:).name})';
        % #########Read Original Image#########
        filesNO(ci) = size(targetPack,1);
        PeriodNum(ci) = floor(filesNO(ci)/ChNum); %Number of imaging periods for each camera
        % Delete unused images (for the case where the red and green channels are not excited simutaneously)
        if filesNO(ci) == ExTotalNum && ~isempty(find(ChAnalyse(ci,:)==0))
            %         if filesNO(ci) >70000 %Manually modify
            tic,
            disp(['Deleting invalid images......']);
            for Chfi = find(ChAnalyse(ci,:)==0) %1:ChNum  %Read the images in the same camera recording folder
                parfor Prdi = 1:PeriodNum(ci)
                    if mod(Prdi,1000) == 0
                        disp(['Deleting: ',num2str(Prdi)]),
                    end
                    FileName = targetsNameTemp{(Prdi-1)*ChNum+Chfi};
                    try
                        delete(FileName);
                    catch
                        disp('File not found');
                    end
                end
            end
            disp(['Deleting finished!']);
            toc,
            % Get the updated list of image files after deleting unused images
            targetPack = dir(['*.',filetype]);
            targetsName(:,ci) = ({targetPack(:).name})';
            filesNO(ci) = size(targetPack,1);
        else
            disp('No images for deleting');
            targetsName(:,ci) = targetsNameTemp;
        end
        PeriodNum(ci) = floor(filesNO(ci)/ChAnlsNum(ci)); %Number of imaging periods for each camera
        PeriodNumRange = PeriodNumRangeList(fi,1):min([PeriodNumRangeList(fi,2),PeriodNum(ci)]); % The actual range of image frames to be processed
        PeriodNumR(ci) = length(PeriodNumRange);
        % Add this section to correct the bug when the duration of Spike2 recording is shorter than the actual imaging duration
        if BlockFrameNum<PeriodNum(ci)
            disp('Bug! The duration of Spike2 recording is shorter than the actual imaging duration');
            ImagingTimeM0 = ImagingTimeM;
            ImagingTimeM(BlockFrameNum+1:PeriodNum(ci)) = ImagingTimeM(end)+max(Period_T)*(1:PeriodNum(ci)-BlockFrameNum);
            disp('Imaging Time Corrected！');
            BlockFrameNum = size(ImagingTimeM,2);
        end
    end
    ImagingTimeM0 = ImagingTimeM;
    ImagingTimeM = ImagingTimeM(PeriodNumRange);
    
    if StimRspFlag % Adjust the values of IdxStartAll and IdxEndAll based on the range of image frames to be processed
        IdxStartAll0 = IdxStartAll;
        IdxEndAll0 = IdxEndAll;
        IdxStartAll = IdxStartAll-PeriodNumRange(1)+1;
        IdxEndAll = IdxEndAll-PeriodNumRange(1)+1;
        BadIdx = find(IdxEndAll>PeriodNumRange(end));
        IdxStartAll(BadIdx) = [];
        IdxEndAll(BadIdx) = [];
    end
    %% Process the results recorded by different cameras
    BlockNum = ceil(PeriodNumR(1)/BlockSz); % Number of groups for grouped data
    IMseriesTave =  zeros(pixR,pixC,CameraNum,ChAnlsNum(1),BlockNum,'uint16');
    IMseriesTmin = IMseriesTave;
    IMseriesTmax = IMseriesTave;
    IMseriesRawAve = zeros(pixR,pixC,CameraNum,ChAnlsNum(1),'uint16');
    IMseriesRawMin = IMseriesRawAve;
    IMseriesAve  = zeros([ceil([pixR*IMscale,pixC*IMscale]),CameraNum,ChAnlsNum(1)],'uint16');
    
    % Generate IMseries
    fdir = dir([SaveFolder1 '\IMseries*.datbin']);
    if ~isempty(fdir)
        % load IMseries
        filename = fullfile(fdir(end).folder,fdir(end).name);
        if strfind(filename,'baseCor.datbin')
            disp('Loading IMseries after baseline correction');
            baseCorFinishFlag = 1;
        else
            disp('Loading IMseries');
            baseCorFinishFlag = 0;
        end
        IMseries = F_load_IMseries_from_datbin_file_Ver1(filename);  % Quickly read from the stored binary file
    else
        baseCorFinishFlag = 0;
        IMseries = zeros([ceil([pixR*IMscale,pixC*IMscale]),PeriodNumR(ci),ChAnlsNum(ci),CameraNum],'uint16');
        shiftsRcd = [];
        tStart = tic;
        for ci = 1:CameraNum
            TempTitle = SessionChTitle{fi,ci};
            %         Idx = find(TempTitle=='_',1,'last');
            TargetPath = TargetList{ci+(fi-1)*CameraNum}; % The folder name of each camera recording is in the format of XXX_G_1, and all the raw data is placed according to the session
            cd(TargetPath);
            disp(['Importing data folder: ',TargetPath]),
            for bi = 1:BlockNum
                disp(['Importing block: ',num2str(bi)]),
                % Read and resize valid images
                BlockFrameS = (bi-1)*BlockSz+1; % Read the starting frame of images in the current block, numbered from 1
                BlockFrameE = min(bi*BlockSz,PeriodNumR(ci)); % Read the ending frame of images in the current block, numbered from 1
                BlockFrameNum = BlockFrameE-BlockFrameS+1; % Total number of frames in the current block
                IMseriesT = zeros([ceil([pixR,pixC]),BlockFrameNum,ChAnlsNum(ci)],'uint16');
                for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %Read the images in the same camera recording folder
                    disp(['Importing(laser): ',ChTag{ci,Chfi}]),
                    IdxCamera = find(CorRatio.Camera==CameraTag{ci}); % Camera index, used to find the corresponding CorRatio.data
                    IdxLaser = find(CorRatio.Laser==ChTag{ci,Chfi}); % Camera index, used to find the corresponding CorRatio.data
                    ImCorRatio = imresize(CorRatio.data(:,:,IdxCamera,IdxLaser),1/ImBinning,'Method','nearest');   % Corresponding camera and laser condition，CorRatio.data
                    tic,
                    for Prdi = BlockFrameS:BlockFrameE
                        if mod(Prdi,1000) == 0
                            disp(['Importing: ',num2str(Prdi)]),
                        end
                        %                 FileName = targetsName{(Prdi-1)*ChNum+ChAnalyseIdx(Chfi),ci};
                        FileName = targetsName{(PeriodNumRange(Prdi)-1)*ChAnlsNum(ci)+Chfi,ci};
                        % Read original image
                        ImgRaw = double(imread(FileName)-RefDark(:,:,ci)); % Read image and substrate camera offset
                        % Correction of illumination uniformity
                        ImgCor = ImgRaw.*ImCorRatio; % uniform/flatted correction
                        %                 figure,
                        %                 PixMin = min(ImgRaw,[],'all');
                        %                 PixMax = max(ImgRaw,[],'all');
                        %                 subplot(1,3,1);imshow(ImgRaw,[PixMin PixMax]),title('Original Image'); %Show original image
                        %                 subplot(1,3,2);imshow(ImCorRatio,[min(ImCorRatio(:)) max(ImCorRatio(:))]),title('Calibration Parameter'); %Display the superpixel segmentation image
                        %                 subplot(1,3,3);imshow(ImgCor,[PixMin PixMax]),title('Image Calibration'); %Display the superpixel segmentation image
                        %                 sgtitle(strrep(FileName,'_','\_'));
                        IMseriesT(:,:,Prdi-BlockFrameS+1,Chfi) = ImgCor;
                    end
                    toc,
                    % Translation correction
                    if MotionCorFlag
                        tic,
                        IMseriesTs = squeeze(IMseriesT(:,:,:,Chfi));
                        disp('Starting motion correction...');
                        if size(shiftsRcd,1) == PeriodNumR(ci) % If there is complete shiftsRcd data, use it directly
                            shifts = shiftsRcd(BlockFrameS:BlockFrameE,:);
                            IMseriesT(:,:,:,Chfi) = apply_shifts(IMseriesTs,shifts,NoRMCorOptions,col_shift);
                        else
                            if Chfi == 1
                                NoRMCorOptions = NoRMCorreSetParms('d1',size(IMseriesTs,1),'d2',size(IMseriesTs,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
                                if ~exist('NoRMCorTemp','var')
                                    [IMseriesT(:,:,:,Chfi),shifts,NoRMCorTemp,NoRMCorOptions,col_shift] = normcorre(IMseriesTs,NoRMCorOptions);
                                else
                                    disp('Correction template alread existed')
                                    [IMseriesT(:,:,:,Chfi),shifts,NoRMCorTemp,NoRMCorOptions,col_shift] = normcorre(IMseriesTs,NoRMCorOptions,NoRMCorTemp);
                                end
                                shiftsRcd = cat(1,shiftsRcd,shifts);
                            else
                                IMseriesT(:,:,:,Chfi) = apply_shifts(IMseriesTs,shifts,NoRMCorOptions,col_shift);
                            end
                        end
                        toc,
                        disp('------------------------------------------------------------------------');
                    end
                    IMseriesTave(:,:,ci,Chfi,bi) = mean(IMseriesT(:,:,:,Chfi),3);
                    IMseriesTmin(:,:,ci,Chfi,bi) = min(IMseriesT(:,:,:,Chfi),[],3);
                    IMseriesTmax(:,:,ci,Chfi,bi) = max(IMseriesT(:,:,:,Chfi),[],3);
                    figure(100),imshow(squeeze(IMseriesTave(:,:,ci,Chfi,bi)),[]);
                    title('Averaged image');
                end
                IMseries(:,:,BlockFrameS:BlockFrameE,:,ci) = imresize(IMseriesT,IMscale,'bilinear');
                %         figure,imshow(mean(IMseries(:,:,BlockFrameS:BlockFrameE,1),3),[]);
                %          figure,imshow(mean(IMseries(:,:,1:1000,1),3),[]);
            end
            clear IMseriesTs IMseriesT;
            % plot shifts
            %             shifts_r = squeeze(cat(3,shiftsRcd(:).shifts));
            %             figure;
            %             ax1 = subplot(211); plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            %             set(gca,'Xtick',[])
            %             ax2 = subplot(212); plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            %             xlabel('timestep','fontsize',14,'fontweight','bold')
            %             linkaxes([ax1,ax2],'x');
            
            
            for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %Read the images in the same camera recording folder
                IMseriesAve(:,:,ci,Chfi) = mean(IMseries(:,:,:,Chfi,ci),3);
                FileName = fullfile(SaveFolder2,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},'_Ave.tif']);
                obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
                WriteIMG(obj,IMseriesAve(:,:,ci,Chfi)'); % save averaged image
                close(obj);
                
                IMseriesRawAve(:,:,ci,Chfi) = mean(IMseriesTave(:,:,ci,Chfi,bi),5);
                FileName = fullfile(SaveFolder2,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},'_RawAve.tif']);% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
                obj = Fast_BigTiff_Write(FileName,pixelsize0,compression);
                WriteIMG(obj,IMseriesRawAve(:,:,ci,Chfi)'); % save averaged image
                close(obj);
                
                IMseriesRawMin(:,:,ci,Chfi) = min(IMseriesTmin(:,:,ci,Chfi,bi),[],5);
                FileName = fullfile(SaveFolder2,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},'_RawMin.tif']);% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
                obj = Fast_BigTiff_Write(FileName,pixelsize0,compression);
                WriteIMG(obj,IMseriesRawMin(:,:,ci,Chfi)'); % save averaged image
                close(obj);
            end
        end
        tEnd = toc(tStart);
        
        %%  registration of different channels
        RegisterMat = cell(CameraNum,1);
        if RegisterFlag
            [optimizer, metric] = imregconfig('multimodal');%The parameter 'modality' specifies the relationship between the fixedIM image and the movingIM image, with two options: 'monomodal' and 'multimodal'
            optimizer.InitialRadius = optimizer.InitialRadius/3.5;%Change the step size of the optimizer to achieve more precise transformation
            % optimizer.MaximumIteCorns = 300;% Change the maximum number of iterations
        end
        
        % Calculate the registration matrix based on the reference image
        if RegisterFlag
            fixedIM = IMseriesAve(:,:,fixCh,1);
            RfixedIM = imref2d(size(fixedIM));
            figure,imshow(fixedIM,[]);
            title(['fixedIM from Camera ',CameraTag{fixCh}]);
            for ci = 1:CameraNum
                if ci~=fixCh
                    movingIM = IMseriesAve(:,:,ci,1);
                    figure,imshow(movingIM,[]);
                    title(['movingIM of Camera ',CameraTag{ci}]);
                    figure(100),imshowpair(movingIM*3,fixedIM,'falsecolor'),title('Before registration');
                    tformSimilarity = imregtform(movingIM,fixedIM,'similarity',optimizer,metric);
                    movingIMReg = imwarp(movingIM,tformSimilarity,'OutputView',RfixedIM);
                    figure(101),imshowpair(movingIMReg,fixedIM,'falsecolor'),title('After registration');
                    %   saveas(gcf,[savefolder1,'\RegistCorn of ',imname(1:12),'.jpg']);
                    RegisterMat{ci,1} = tformSimilarity;
                    RefDarkRegister(:,:,ci) = imwarp(RefDark(:,:,ci),tformSimilarity,'OutputView',RfixedIM);
                    
                    tic
                    IMseries(:,:,:,:,ci) = imwarp(IMseries(:,:,:,:,ci),tformSimilarity,'OutputView',RfixedIM);
                    toc,
                    for Chfi = 1:ChAnlsNum(ci)
                        IMseriesAve(:,:,ci,Chfi) = mean(IMseries(:,:,:,Chfi,ci),3);
                        FileName = fullfile(SaveFolder2,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},'_AveRgst','.tif']);% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
                        obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
                        WriteIMG(obj,IMseriesAve(:,:,ci,Chfi)'); % save averaged image
                        close(obj);
                    end
                    disp('Multi camera channel registeration finished!');
                end
            end
        end
        
        %% Save IMseries
        if forceSaveFlag
            tic,
            SzIMseries = size(IMseries);
            if length(SzIMseries)==4
                SzIMseries = [SzIMseries,1];
            end
            fn = sprintf('IMseries1_%d_%d_%d_%d_%d_uint16_%s.datbin',SzIMseries,SessionTitle{fi,1});
            disp(['saving ',fn,'......']);
            FileName = fullfile(SaveFolder1,fn);
            f = fopen(FileName, 'w');
            fwrite(f, IMseries, 'uint16');
            fclose(f);
            disp('IMseries saved');
            toc;
        end
    end
    [ImSz1,ImSz2,ImSz3,ImSz4,ImSz5] = size(IMseries);
    %% #################### Load WholeCortexMask and ROIs in ImageJ ####################
    % load WholeCortexMask
    cd(SaveFolder);
    %     MaskName = '\WholeCortexMask.tif';
    MaskName = [SaveFolder,'\WholeCortexMask.tif'];
    if ~exist(MaskName,'file')
        uiwait(msgbox('Please generate mask and ROIs','Warning','modal'));
    end
    WholeCortexMask0 = imread(MaskName);
    WholeCortexMask0(find(WholeCortexMask0 ==  0)) = NaN;
    WholeCortexMask0(find(WholeCortexMask0 == 255)) = 1;
    WholeCortexMask = imresize(WholeCortexMask0,IMscale,'Method','bilinear');
    %          WholeCortexMask = imresize(WholeCortexMask0,IMscale,'Method','bicubic');
    %          WholeCortexMask = imresize(WholeCortexMask0,IMscale,'Method','bilinear');
    figure,imshow(WholeCortexMask,[]);
    title('Mask imresized');
    %     [ImSz1m,ImSz2m] = size(WholeCortexMask);
    for pixCi = 2:ImSz1-1
        for pixRi = 2:ImSz2-1
            nhoods = WholeCortexMask(pixCi-1:pixCi+1,pixRi-1:pixRi+1);
            if length(find(nhoods== WholeCortexMask(pixCi,pixRi)))==1
                WholeCortexMask(pixCi,pixRi) = nhoods(1,1);
            end
        end
    end
    figure,imshow(WholeCortexMask,[]);
    title('Mask imresized and filled');
    FileName = [MaskName(1:end-4),'_reScale.tif'];
    obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
    WriteIMG(obj,WholeCortexMask'); % save averaged image
    close(obj);
    
    % load ROIs
    ROI_Name = [SaveFolder,'\RoiSet.zip'];
    %     ROIs = poly2mask(ROIs.mnCoordinates(:, 1), ROIs.mnCoordinates(:, 2),ImSz1,ImSz2);
    ROIs = ReadImageJROI(ROI_Name);
    for ri = 1:size(ROIs,2)
        ROI = ROIs{ri};
        switch ROI.strType
            case 'Rectangle'
                ROIxy = ROI.vnRectBounds+1;
                mask0 = uint8(zeros(pixR,pixC));
                %                 mask0 = logical(zeros(pixR,pixC));
                mask0(ROIxy(1):ROIxy(3),ROIxy(2):ROIxy(4)) = 1;
            case 'Polygon'
                ROIxy = ROI.mnCoordinates+1;
                mask0 = uint8(poly2mask(ROIxy(:, 1), ROIxy(:, 2),pixR,pixC));
            otherwise
                disp('ROI type ERROR!');
        end
        mask = imresize(mask0,IMscale,'Method','nearest');
        if ri>1
            mask = mask.*WholeCortexMask;
            mask0 = mask0.*WholeCortexMask0;
        end
        [mask_r{ri}, mask_c{ri}] = find(mask > 0);
        mask_id{ri} = find(mask > 0);
        figure,imshow(mask0,[]);
        title(strrep(ROI.strName,'_','\_'));
        close;
        figure,imshow(mask,[]);
        title(['Resized: ',strrep(ROI.strName,'_','\_')]);
        close;
    end
    %% Frame number of the basal image
    switch BasalFlagOfSession % Methods for generating the baseline image for the whole session: 1-Specified frame；2-Based on the stimulus marker；3. The X% lower value for each pixel in the entire session; 4-The union set of baseline from the begin of session and each trial; 5-based on other labelled baseline, such as REM during sleep wake cycles
        case 1 % 1-Specified frame
            IdxBaseStart = IdxBaseList(fi,1); % Starting frame of the entire session's baseline
            IdxBaseEnd = IdxBaseList(fi,2); % Ending frame of the entire session's baseline
            if ~isnan(IdxBaseStart) & ~isnan(IdxBaseEnd) % If the baseline is not specified, use the previous IMseriesBasalIM
                IdxStart = max([IdxBaseStart,PeriodNumRange(1)])-PeriodNumRange(1)+1;
                IdxEnd = min([IdxBaseEnd,PeriodNumRange(end)])-PeriodNumRange(1)+1;
                if IdxBaseEnd==0
                    IdxEnd = ImSz3;
                end
            else
                IdxStart = IdxStart_bkp{fi-1}; % Using the baseline image from previous session
                IdxEnd = IdxEnd_bkp{fi-1}; % Using the baseline image from previous session
            end
            
        case 2  % 2-Based on the stimulus marker. The time before the first stimulus is defined as the end time of the baseline for the entire session.
            Baseline0Start = 0; % The starting time of the baseline for the entire session.
            Baseline0End = Stim.times(1,1); % Time of the first stimulus
            IdxStart = ceil(find(Imaging.times>=Baseline0Start,1)/ChNum/2);
            IdxEnd = floor(find(Imaging.times<=Baseline0End,1,'last')/ChNum/2);
            
        case 3  % 3-The X% lower value for each pixel in the entire session
            IdxBaseStart = IdxBaseList(fi,1); % Starting frame of the entire session's baseline
            IdxBaseEnd = IdxBaseList(fi,2); % Ending frame of the entire session's baseline
            IdxStart = max([IdxBaseStart,PeriodNumRange(1)])-PeriodNumRange(1)+1;
            IdxEnd = min([IdxBaseEnd,PeriodNumRange(end)])-PeriodNumRange(1)+1;
            
        case 4 % 4-The union set of baseline from the begin of session and each trial
            IdxStart = cat(1,IdxBaseList(fi,1),IdxStartAll);
            IdxEnd = cat(1,IdxStartAll(1)-IdxCorLen,IdxStimAll-IdxCorLen);
            
        case 5 % 5-based on other labelled baseline, such as REM during sleep wake cycles
            % Analysis by AccuSleep (1 = REM sleep, 2 = wakefulness, 3 = NREM sleep, 4 = undefined)
            StimRcdPath1 = fullfile(StimRcdPath,'SleepState');
            FileName = fullfile(StimRcdPath1,'SleepStateLabel.mat');
            load(FileName);
            
            %     labelsLen = size(labels,1);
            %     labelsTimeM = size(labels,1)*EpochTime;
            
            ImStamp.Flag = zeros(PeriodNumR(ci),1)*nan;
            ImStamp.Note = ImStamp.Flag;
            StimModeCount = zeros(max(labels),1);
            for Prdi = 1:PeriodNumR(ci) % Process each frame
                try
                    StimMode = labels(ceil(ImagingTimeM(Prdi)/EpochTime)); % If the last frame does not have a corresponding sleep state recording, skip the error and use the previous state.
                    ImStamp.Flag(Prdi,1) = StimMode;
                    StimModeCount(StimMode,1) = StimModeCount(StimMode,1)+1;
                    ImStamp.Note(Prdi,1) = StimModeCount(StimMode,1); % The number of specific trial that belongs to the event or stimulus.
                    ImStamp.Note(Prdi,2) = Prdi; % The total number of trials.
                end
                IdxStart = NaN;
                IdxEnd = NaN;
                FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_ImStamp.mat']);% Save parameters of the linear fit : k (slope) and b (intercept).
                save(FileName,'ImStamp');
            end
            
            if BasalFlagOfSession == 5
                IdxBaseFrame = find(ImStamp.Flag==1); % For data during sleep-wake cycles, use the REM state as the baseline.
            else
                IdxBaseFrame = [];
                for j = 1:size(IdxStart,1)
                    IdxBaseFrame = cat(2,IdxBaseFrame,IdxStart(j):IdxEnd(j));
                end
            end
            IdxBaseFrame_bkp{fi} = IdxBaseFrame;
            IdxStart_bkp{fi} = IdxStart;
            IdxEnd_bkp{fi} = IdxEnd;
            if baseCorFinishFlag
                disp('Background has been already subtracted!');
                disp('Baseline has been already corrected!');
            else
                %% Substract background
                BackgroundRcd = zeros(CameraNum,ChAnlsNum(1),PeriodNumR(1),'uint16');
                ChNo = 0;
                figure,title('BackgroundRcd');hold on;
                for ci = 1:CameraNum
                    for Chfi = 1:ChAnlsNum(ci) %  Different excitation lights for the same camera
                        IMseriesT = IMseries(:,:,:,Chfi,ci);
                        IMseriesT = reshape(IMseriesT,[],PeriodNumR(ci));
                        %     BackgroundRcd(ci,Chfi,:) = uint16(min(IMseriesT(mask_id{1},:),[],1)*BgToAutoF(ci,Chfi));
                        BackgroundRcd(ci,Chfi,:) = uint16(mean(IMseriesT(mask_id{1},:),1));
                        IMseriesT = IMseries(:,:,:,Chfi,ci)-BackgroundRcd(ci,Chfi,:);
                        IMseries(:,:,:,Chfi,ci) = IMseriesT.*uint16(WholeCortexMask);  % Set the values in pixels of the vascular region to 0
                        %             IMseriesAve(:,:,ci,Chfi) = mean(IMseries(:,:,:,Chfi,ci),3);
                        %             FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},'_AveSubBG','.tif']);
                        %             obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
                        %             WriteIMG(obj,IMseriesAve(:,:,ci,Chfi)'); % save averaged image
                        %             close(obj);
                        plot(squeeze(BackgroundRcd(ci,Chfi,:)));
                        ChNo = ChNo+1;
                        legendTag{ChNo,1} = [CameraTag{ci},'_',ChTag{ci,Chfi}];
                    end
                end
                legend(strrep(legendTag,'_','\_'));
                disp('Background subtracted!');
                %% baseline correction (such as for photobleaching)
                if BaselineCorFlag
                    x = 1:ImSz3;
                    BaseRange = 1:round(length(x)*0.05); % Use the baseline fit obtained from the first 5% time period as the expected true value for the baseline.
                    fdir = dir([MfileDir '\Funs\createFits.m']);
                    if ~isempty(fdir)
                        disp('createFits.m file exists');
                    else
                        %         yVarNameList = cell(CameraNum,ChAnlsNum(1));
                        if BaselineCorFlag == 1
                            for ci = 1:CameraNum
                                for Chfi = 1:ChAnlsNum(ci)
                                    IMseriesT = IMseries(:,:,:,Chfi,ci);
                                    IMseriesT = reshape(IMseriesT,[],size(IMseriesT,3));
                                    y = mean(IMseriesT(mask_id{2},:),1);
                                    yVarName = ['y',CameraTag{ci},ChTag{ci,Chfi}];
                                    %                 yVarNameList{ci,Chfi} = yVarName;
                                    eval([yVarName,'=y;']);
                                    eval(['cftool(x,',yVarName,');']);
                                end
                            end
                            % Please manually save the automatically generated createFits.m file from MATLAB to the mfile\Funs directory.
                            helpdlg('Please generate createFits.m!');
                            pause;
                            uiwait(msgbox({'Please use sftool generate a function and saveas "createFits"';'Make sure ';'Finished?'},'Warning','modal'));
                        end
                    end
                    
                    % baseline correction for each pixel
                    tic,
                    disp('Baseline correction starting...');
                    BlockParNum = 10; % Number of parallel computing groups (grouping to reduce the memory consumption of each parfor iteration).
                    BlockParSz = ceil(ImSz2/BlockParNum);
                    PixiList = find(any(WholeCortexMask,2));
                    PixiList = PixiList';
                    for Pixi = PixiList
                        disp(['pixel row: ',num2str(Pixi),'/',num2str(ImSz1)]);
                        for Pbi = 1:BlockParNum
                            PixjList = (Pbi-1)*BlockParSz+1:min(Pbi*BlockParSz,ImSz2);
                            IMseriesT = squeeze(IMseries(Pixi,PixjList,:,:,:));
                            IMseriesTcor = IMseriesT;
                            parfor PixjIdx = 1:length(PixjList)
                                %                             Pixi = 121; % for debug
                                %                             Pixj = 210;
                                %                             PixjIdx = 18;
                                Pixj = PixjList(PixjIdx);
                                if WholeCortexMask(Pixi,Pixj)
                                    %                             tic,
                                    disp(['Pix row:',num2str(Pixi),', Pix column:',num2str(Pixj)]);
                                    y = zeros(4,PeriodNumR(1));
                                    for ci = 1:CameraNum
                                        for Chfi = 1:ChAnlsNum(ci)
                                            Idx = (ci-1)*CameraNum+Chfi;
                                            y(Idx,:) = double(squeeze(IMseriesT(PixjIdx,:,Chfi,ci)));
                                        end
                                    end
                                    if BaselineCorFlag == 1
                                        [fitresult, gof] = createFits(x,  y(1,:),  y(2,:),  y(3,:), y(4,:));
                                    end
                                    yCorPix = zeros(PeriodNumR(1),ChAnlsNum(ci),CameraNum);
                                    for ci = 1:CameraNum
                                        for Chfi = 1:ChAnlsNum(ci)
                                            Idx = (ci-1)*CameraNum+Chfi;
                                            switch BaselineCorFlag
                                                case 1
                                                    FitC = coeffvalues(fitresult{Idx});
                                                    yFit = FitC(1)*exp(FitC(2)*x) + FitC(3)*exp(FitC(4)*x);
                                                case 2
                                                    [BEADSx, yFit, BEADScost] = F_beads_Ver0(y(Idx,:),0);
                                                case 3
                                                    yFit = smooth(y(Idx,:),0.05,'rloess');
                                                    yFit = yFit';
                                            end
                                            yCor = uint16(y(Idx,:)./yFit*mean(yFit(BaseRange)));
                                            yCorPix(:,Chfi,ci) = yCor;
                                            %                                     figure,plot(1:PeriodNumR(1),y(Idx,:),1:PeriodNumR(1),yFit,1:PeriodNumR(1),yCor);
                                            %                                     legend('y','yFit','yCor');
                                            %                                     title(['Trace of Pix row:',num2str(Pixi),', Pix column:',num2str(Pixj),', Camera:',num2str(ci),', Ex: ',num2str(Chfi)]);
                                        end
                                    end
                                    IMseriesTcor(PixjIdx,:,:,:) = yCorPix;
                                    y = [];
                                    fitresult = [];
                                    yFit = [];
                                    yCor = [];
                                    yCorPix = [];
                                    %                             toc;
                                end
                            end
                            IMseries(Pixi,PixjList,:,:,:) = IMseriesTcor;
                            clear IMseriesT IMseriesTcor;
                        end
                        % Clear the parallel pool when the RAM usage is too high
                        if mod(Pixi,10)==0
                            [userview,systemview] = memory;
                            RAMratio = systemview.PhysicalMemory.Available/systemview.PhysicalMemory.Total;
                            if RAMratio<0.3
                                disp('Clean parfor pool');
                                delete(gcp('nocreate'));
                            end
                        end
                    end
                    toc
                    
                    % Save IMseries
                    if forceSaveFlag
                        tic,
                        SzIMseries = size(IMseries);
                        if length(SzIMseries)==4
                            SzIMseries = [SzIMseries,1];
                        end
                        fn = sprintf('IMseries1_%d_%d_%d_%d_%d_uint16_%s_baseCor.datbin',SzIMseries,SessionTitle{fi,1});
                        disp(['saving ',fn,'......']);
                        FileName = fullfile(SaveFolder1,fn);
                        f = fopen(FileName, 'w');
                        fwrite(f, IMseries, 'uint16');
                        fclose(f);
                        disp('IMseries saved');
                        toc;
                    end
                end
            end
            % Temporarily generate basal images for subsequent calculations
            IMseriesBasalIM = zeros(ImSz1,ImSz2,ChAnlsNum(1),CameraNum,'uint16');
            LowThd = 0.1; % The threshold for filtering whether each pixel point can be used for baseline calculation
            for ci = 1:CameraNum
                for Chfi = 1:ChAnlsNum(ci)
                    [IMseriesBasalIM(:,:,Chfi,ci),IMseriesBasalIM_Thd(:,:,Chfi,ci)] = F_BasalIM_Ver1(IMseries(:,:,:,Chfi,ci),BasalFlagOfSession,LowThd,IdxBaseFrame);
                    %         figure,imshow(IMseriesBasalIM(:,:,Chfi,ci),[]);
                    % title(['Basal image: ',num2str(ci),', ',num2str(Chfi)]);
                end
            end
            %% Spectra unmixing
            if SpectraUnmixFlag & CameraNum>1
                tic,
                disp('Spectra Unmixing starting...');
                for Chfi = 1:ChAnlsNum(ci)
                    %                         img = mean(IMseries(:,:,:,Chfi,1),3);
                    %                         figure,imshow(img,[]);
                    %                         title('Before Spectra Unmixing');
                    
                    %             measuredR = IMseries(:,:,:,Chfi,2)-IMseriesBasalIM(:,:,Chfi,2);
                    dG = single(IMseries(:,:,:,Chfi,1))-single(IMseriesBasalIM(:,:,Chfi,1));
                    dR = single(IMseries(:,:,:,Chfi,2))-single(IMseriesBasalIM(:,:,Chfi,2));
                    IMseries(:,:,:,Chfi,1) = uint16((dG - Crg*dR)/(1-Cgr*Crg)+single(IMseriesBasalIM(:,:,Chfi,1)));
                    IMseries(:,:,:,Chfi,2) = uint16((dR - Cgr*dG)/(1-Cgr*Crg)+single(IMseriesBasalIM(:,:,Chfi,2)));
                    %                         img = mean(IMseries(:,:,:,Chfi,1),3);
                    %                         figure,imshow(img,[]);
                    %                         title('After Spectra Unmixing');
                end
                clear dG dR;
                disp('Spectra Unmixing Finished!');
                toc,
            end
            
            
            % Generate basal image and AveSubBG image
            IMseriesBasalIM = zeros(ImSz1,ImSz2,ChAnlsNum(1),CameraNum,'uint16');
            for ci = 1:CameraNum
                for Chfi = 1:ChAnlsNum(ci)
                    [IMseriesBasalIM(:,:,Chfi,ci),IMseriesBasalIM_Thd(:,:,Chfi,ci)] = F_BasalIM_Ver1(IMseries(:,:,:,Chfi,ci),BasalFlagOfSession,LowThd,IdxBaseFrame);
                    %         figure,imshow(IMseriesBasalIM(:,:,Chfi,ci),[]);
                    % title(['Basal image: ',num2str(ci),', ',num2str(Chfi)]);
                    FileName = fullfile(SaveFolder2,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},'_Basal','.tif']);
                    obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
                    WriteIMG(obj,IMseriesBasalIM(:,:,Chfi,ci)'); % save averaged basal image
                    close(obj);
                    
                    IMseriesAve(:,:,ci,Chfi) = mean(IMseries(:,:,:,Chfi,ci),3);
                    FileName = fullfile(SaveFolder2,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},'_AveSubBG','.tif']);
                    obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
                    WriteIMG(obj,IMseriesAve(:,:,ci,Chfi)'); % save averaged image after subtract background
                    close(obj);
                end
            end
            IMseriesBasalIM_bkp{fi} = IMseriesBasalIM;
            %% Calculate response
            for ci = 1:CameraNum
                %% Artifacts correction
                tic,
                switch CorrectMethodFlag
                    case 1
                        % GB or RB Cor images
                        disp('Artifacts correction by ratio image...');
                        Chfi = 1; % Excited by 488 or 561 nm
                        ChNormIdx = 2; % Excited by 405 nm
                        %             IMseriesMCor(:,:,:,Chfi,:) = double(IMseriesM(:,:,:,Chfi,:))./double(IMseriesM(:,:,:,ChNormIdx,:));
                        IMseriesCor = single(IMseries(:,:,:,Chfi,ci))./single(IMseries(:,:,:,ChNormIdx,ci));
                    case 2
                        SaveFolder3 = [SaveFolder1,'\Regression(405nm) evaluation'];
                        mkdir(SaveFolder3);
                        disp('Artifacts correction by fitted 405 channel...');
                        if BasalFlagOfSession == 5
                            IdxFit = find(ImStamp.Flag==2); % 对睡眠觉醒数据，以Wake对应时刻为fit的数据点
                        else
                            IdxFit = IdxBaseFrame; % Correction based on the baseline image
                        end
                        
                        IMseriesT = single(squeeze(IMseries(:,:,IdxFit,:,ci)));
                        IMseriesT(find(IMseriesT==0)) = NaN;
                        %                     figure,imshow(mean(IMseriesT(:,:,:,1),3),[]);
                        %                 figure,imshow(mean(IMseriesT(:,:,:,2),3),[]);
                        FilterSigma = 2;
                        FilterWidth = 2*ceil(2*FilterSigma)+1;
                        ImFilter=fspecial('gaussian',FilterWidth,FilterSigma);
                        for i = 1:length(IdxFit)
                            %                     disp(i);
                            img = squeeze(IMseriesT(:,:,i,1));
                            IMseriesT(:,:,i,1) = nanconv(img,ImFilter, 'nanout');
                            img = squeeze(IMseriesT(:,:,i,2));
                            IMseriesT(:,:,i,2) = nanconv(img,ImFilter, 'nanout');
                        end
                        %                 figure,imshow(mean(IMseriesT(:,:,:,1),3,'omitnan'),[]);
                        %                 figure,imshow(mean(IMseriesT(:,:,:,2),3,'omitnan'),[]);
                        thetaIM = zeros(ImSz1,ImSz2,2);
                        for Pixi = 1:ImSz1
                            disp(['pixel row: ',num2str(Pixi),'/',num2str(ImSz1)]);
                            for Pixj = 1:ImSz2
                                %                         Pixi = 200;
                                %                         Pixj = 200;
                                if WholeCortexMask(Pixi,Pixj)
                                    PixV = squeeze(IMseriesT(Pixi,Pixj,:,:));
                                    %                             PixV = squeeze(IMseries(Pixi,Pixj,IdxFit,:,ci));
                                    PixV = double(PixV');
                                    PixVB = PixV(1,:)';
                                    PixVA = [PixV(2,:)' ones(length(IdxFit),1)];
                                    theta = PixVA\PixVB;
                                    thetaIM(Pixi,Pixj,1:2) = theta;
                                    if mod(Pixi,50)==0 & mod(Pixj,50)==0 % random inspection of pixels to evaluate the fit performance
                                        % Visualize for debug
                                        InputTag = ['pixel r',num2str(Pixi),'c',num2str(Pixj)];
                                        PixVall = zeros(ImSz3,2);
                                        PixVall = double(squeeze(IMseries(Pixi,Pixj,:,:,ci))');
                                        PixVAfit = theta(1)*PixVall(2,:)+theta(2);
                                        figure(101),clf,hold on;
                                        %                             plot(ImagingTimeM(IdxFit),PixVall(1,:),'g'); plot(ImagingTimeM(IdxFit),PixVall(2,:),'b'); plot(ImagingTimeM(IdxFit),PixVAfit,'k');
                                        plot(ImagingTimeM,PixVall(1,:),'g'); plot(ImagingTimeM,PixVall(2,:),'b'); plot(ImagingTimeM,PixVAfit,'k');
                                        set(gca,'FontSize',10); xlim([0 ImagingTimeM(end)]);
                                        legend('488 nm','405 nm','Fitted 405 nm');
                                        xlabel('Time (sec)','FontSize',12);
                                        ylabel('Fluorescence (au)','FontSize',12);
                                        title('488 nm and "fitted" 405 nm signals','FontSize',12);
                                        FileName =  fullfile(SaveFolder3,[CameraTag{ci},'_',InputTag,'_All signals','.fig']);
                                        saveas(gcf,FileName);
                                        FileName = fullfile(SaveFolder3,[CameraTag{ci},'_',InputTag,'_All signals','.png']);
                                        saveas(gcf,FileName);
                                        FileName = fullfile(SaveFolder3,[CameraTag{ci},'_',InputTag,'_fit cor_pix',num2str(Pixi),' pix',num2str(Pixj),'.mat']);
                                        save(FileName,'ImagingTimeM','PixV','PixVAfit','theta');
                                        
                                        figure(102),clf,hold on;
                                        PixRsp = PixVall(1,:)/mean(PixVall(1,IdxFit))-1;
                                        plot(PixRsp,'g');
                                        PixRspCor = (PixVall(1,:)-PixVAfit)/mean(PixVAfit(IdxFit));
                                        plot(PixRspCor,'k');
                                        %                                         plot(ImagingTimeM(IdxFit),PixVall(1,:),'g'); plot(ImagingTimeM(IdxFit),PixVall(2,:),'b'); plot(ImagingTimeM(IdxFit),PixVAfit,'k');
                                        set(gca,'FontSize',10); xlim([0 length(PixVall)]);
                                        legend('488 nm','Fitted 405 nm');
                                        xlabel('Time (sec)','FontSize',12);
                                        ylabel('dF/F0','FontSize',12);
                                        %                                 ylim([-0.4 0.4]);
                                        title(['dF/F0 signals of ',InputTag],'FontSize',12);
                                        FileName = fullfile(SaveFolder3,[CameraTag{ci},'_',InputTag,'_response','.fig']);
                                        saveas(gcf,FileName);
                                        FileName = fullfile(SaveFolder3,[CameraTag{ci},'_',InputTag,'_response','.png']);
                                        saveas(gcf,FileName);
                                        
                                        figure(103),clf,scatter(PixV(2,:),PixV(1,:),10);hold on;
                                        Xmin = min(PixV(2,:));
                                        Xmax = max(PixV(2,:));
                                        Ymin = min(PixV(1,:));
                                        Ymax = max(PixV(1,:));
                                        % Ymin = theta(1)*Xmin+theta(2);
                                        % Ymax = theta(1)*Xmax+theta(2);
                                        plot([Xmin Xmax],[Ymin Ymax]);
                                        txtPx = (Xmax-Xmin)*0.001+Xmin;
                                        txtPy = Ymax-(Ymax-Ymin)*0.05;
                                        txtT = ['y = ',num2str(theta(1)),'x+',num2str(theta(2))];
                                        text(txtPx,txtPy,txtT,'fontsize',20);
                                        xlabel('MFI(Ex405nm)','FontSize',12);
                                        ylabel('MFI(Ex488nm)','FontSize',12);
                                        title(['Fit of ',InputTag],'FontSize',12);
                                        %                             xlim([1900 2300]);
                                        %                             set(gca,'xtick',1900:100:2300);
                                        %                             ylim([2900 3400]);
                                        %                             set(gca,'ytick',2900:100:3400);
                                        FileName = fullfile(SaveFolder3,[CameraTag{ci},'_',InputTag,'_Fit regression','.fig']);
                                        saveas(gcf,FileName);
                                        FileName = fullfile(SaveFolder3,[CameraTag{ci},'_',InputTag,'_Fit regression','.png']);
                                        saveas(gcf,FileName);
                                        close(101:103)
                                    end
                                end
                            end
                        end
                        thetaIMk = thetaIM(:,:,1);
                        thetaIMb = thetaIM(:,:,2);
                        FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_fit value','.mat']);% 保存线性fit参数k，斜率，参数b，截距
                        save(FileName,'thetaIM');
                        yMin = min(thetaIMk(:));
                        yMax = max(thetaIMk(:));
                        figure,imshow(thetaIMk,[yMin yMax]);
                        title('k value');
                        colormap jet;
                        colorbarNote = colorbar;
                        %                 colorbarNote.Label.String = 'a.u.';
                        %                 colorbarNote.Label.Rotation = 270;
                        %                 colorbarNote.Label.FontSize = 16;
                        colorbarNote.Ticks = [yMin yMax];
                        FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_fit k','.png']);% 线性fit参数k，斜率
                        saveas(gcf,FileName);
                        
                        yMin = min(thetaIMb(:));
                        yMax = max(thetaIMb(:));
                        figure,imshow(thetaIMb,[yMin yMax]);
                        title('b value');
                        colormap jet;
                        colorbarNote = colorbar;
                        colorbarNote.Ticks = [yMin yMax];
                        FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_fit b','.png']);% 线性fit参数b，截距
                        saveas(gcf,FileName);
                        Chfi = 1; % Excited by 488 or 561 nm
                        ChNormIdx =2; % Excited by 405 nm
                        IMseriesCor = single(IMseries(:,:,:,Chfi,ci))-thetaIMk.*single(IMseries(:,:,:,ChNormIdx,ci))-thetaIMb;
                        IMseriesChNormBasalIM = IMseriesBasalIM(:,:,ChNormIdx,ci);
                        IMseriesChNormBasalIMfit = thetaIMk.*single(IMseriesChNormBasalIM)+thetaIMb;
                        figure,
                        subplot(1,2,1),imshow(IMseriesChNormBasalIM,[]);title('Ex 405 nm');
                        subplot(1,2,2),imshow(IMseriesChNormBasalIMfit,[]);title('Ex 405 nm fitted');
                        FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},'_Ex405&fit images','.fig']);% 线性fit参数b，截距
                        saveas(gcf,FileName);
                        IMseriesCor = IMseriesCor+IMseriesChNormBasalIMfit; % To prevent F0 from being too small, add the fitted 405 value after the fit
                        clear IMseriesT;
                end
                toc,
                %% Analysis of response across the whole session
                if WholeSessionRspFlag
                    disp('Analyzing response of the whole session...'),
                    disp(['Generate basal images by method ',num2str(BasalFlagOfSession)]);
                    % Export the basal image after correction
                    Chfi = 1;
                    ChNorm = 2;
                    [IMseriesCorBasalIM(:,:,ci),IMseriesCorBasalIM_Thd(:,:,ci)] = F_BasalIM_Ver1(IMseriesCor,BasalFlagOfSession,LowThd,IdxBaseFrame);
                    IMseriesCorBasalIM_bkp{fi} = IMseriesCorBasalIM(:,:,ci);
                    figure,imshow(IMseriesCorBasalIM(:,:,ci),[]);
                    title('IMseriesCor baseline average');
                    FileName = fullfile(SaveFolder2,[SessionChTitle{fi,1},'_',ChTag{ci,Chfi},ChTag{ci,ChNorm},'cor_Basal','.tif']);
                    obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
                    WriteIMG(obj,IMseriesCorBasalIM(:,:,ci)');
                    close(obj);
                    if exist('IMseriesCorBasalIM_Thd')
                        FileName = fullfile(SaveFolder2,[SessionChTitle{fi,1},'_',ChTag{ci,Chfi},ChTag{ci,ChNorm},'cor_BasalThd','.tif']);% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
                        obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
                        WriteIMG(obj,IMseriesCorBasalIM_Thd(:,:,ci)'); % save averaged image
                        close(obj);
                    end
                    IMseriesCorAve(:,:,ci) = mean(IMseriesCor,3);
                    FileName = fullfile(SaveFolder2,[SessionChTitle{fi,1},'_',ChTag{ci,Chfi},ChTag{ci,ChNorm},'cor_Ave','.tif']);% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
                    obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
                    WriteIMG(obj,IMseriesCorAve(:,:,ci)'); % save averaged image
                    close(obj);
                    
                    % ######### Export the response image for a single channel #########
                    compression = 0;
                    IMseriesRsp = zeros(ImSz1,ImSz2,ImSz3,ImSz4,'single');
                    tic,
                    for Chfi = 1:ChAnlsNum(ci)
                        IMseriesRsp(:,:,:,Chfi) = squeeze(single(IMseries(:,:,:,Chfi,ci)).*single(WholeCortexMask)./single(IMseriesBasalIM(:,:,Chfi,ci)))-1;
                    end
                    IdxBad = find(abs(IMseriesRsp)==Inf);
                    if IdxBad
                        IMseriesRsp(IdxBad) = 0;
                    end
                    for Chfi = 1:ChAnlsNum(ci)
                        if SaveSingleChFlag(ci)
                            FileName = fullfile(SaveFolder1,[SessionChTitle{fi,ci},'_',ChTag{ci,Chfi},'_Rsp.tif']);
                            F_WriteBigTiff_Ver1(FileName,pixelsize,compression,IMseriesRsp(:,:,:,Chfi));
                            if BasalFlagOfSession == 5
                                % Output the average images during different sleep states and .mat files for sleep states
                                StatusStamp = ImStamp.Flag;
                                StatusLabels  = {'REM','Wake','NREM'};
                                IMseries3D = IMseriesRsp(:,:,:,Chfi); %  IMseries3D = IMseriesT;
                                FileNameSuffix = fullfile(SaveFolder1,[SessionChTitle{fi,ci},'_',ChTag{ci,Chfi},'_']);% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
                                IMofStatus = F_StatusMeanIM_Ver1(IMseries3D,StatusStamp,StatusLabels,FileNameSuffix,pixelsize,compression);
                            end
                        else
                            disp('No saving response data based on singel channel (Whole session)');
                        end
                    end
                    toc,
                    % ######### Export the response image after correction #########
                    tic,
                    IMseriesRspCor = IMseriesCor.*single(WholeCortexMask)./IMseriesCorBasalIM(:,:,ci)-1;
                    IdxBad = find(abs(IMseriesRspCor)==Inf);
                    if IdxBad
                        IMseriesRsp(IMseriesRspCor) = 0;
                    end
                    Chfi = 1;
                    FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},ChTag{ci,ChNormIdx},'cor_Rsp','.tif']); % num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
                    F_WriteBigTiff_Ver1(FileName,pixelsize,compression,IMseriesRspCor);
                    if BasalFlagOfSession == 5
                        % Output the average images during different sleep states and .mat files for sleep states
                        StatusStamp = ImStamp.Flag;
                        StatusLabels  = {'REM','Wake','NREM'};
                        FileNameSuffix = fullfile(SaveFolder1,[SessionChTitle{fi,ci},'_',ChTag{ci,Chfi},ChTag{ci,ChNormIdx},'cor_']);% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
                        IMofStatus = F_StatusMeanIM_Ver1(IMseriesRspCor,StatusStamp,StatusLabels,FileNameSuffix,pixelsize,compression);
                    end
                    toc,
                    
                    % response of each ROIs across the whole session
                    % calculate response
                    tic,
                    i = 1;
                    RspOfSession.ChS = cell(i,size(ROIs,2)-1); % Response, excited by 488 or 561 nm
                    RspOfSession.ChNorm = cell(i,size(ROIs,2)-1); % Response, excited by 405 nm
                    RspOfSession.Cor = cell(i,size(ROIs,2)-1); % Response, after correction
                    RspOfSession.StartTimeInSec = ImagingTimeM(1);
                    RspOfSession.TimeInSec = ImagingTimeM'-RspOfSession.StartTimeInSec;
                    RspOfSession.TimeInMin = RspOfSession.TimeInSec/60;
                    RspOfSession.TimeInHrs = RspOfSession.TimeInSec/3600;
                    MFIOfSession.ChS = cell(i,size(ROIs,2)-1); % Average fluorescence intensities, excited by 488 or 561 nm
                    MFIOfSession.ChNorm = cell(i,size(ROIs,2)-1); % Average fluorescence intensities, excited by 405 nm
                    MFIOfSession.Cor = cell(i,size(ROIs,2)-1); % Average fluorescence intensities, after correction
                    
                    % ######## Excited by 488 or 561 nm ########
                    IMseriesRspChS = IMseriesRsp(:,:,:,1);
                    IMseriesRspChS = reshape(IMseriesRspChS,[],PeriodNumR(ci));
                    %             IMseriesRspChS = IMseriesRspChS';
                    % ######## Excited by 405 nm ########
                    IMseriesRspChNorm = IMseriesRsp(:,:,:,2);
                    IMseriesRspChNorm = reshape(IMseriesRspChNorm,[],PeriodNumR(ci));
                    %             IMseriesRspChNorm = IMseriesRspChNorm';
                    % ######## After correction ########
                    IMseriesRspCor = reshape(IMseriesRspCor,[],PeriodNumR(ci));
                    %             IMseriesRspCor = IMseriesRspCor';
                    if MFIrecordFlag  % Record the average fluorescence intensity of the ROI over time
                        % ######## Excited by 488 or 561 nm ########
                        IMseriesMFIChS = IMseries(:,:,:,1,ci);
                        IMseriesMFIChS = reshape(IMseriesMFIChS,[],PeriodNumR(ci));
                        %                 IMseriesMFIChS = IMseriesMFIChS';
                        % ######## Excited by 405 nm ########
                        IMseriesMFIChNorm = IMseries(:,:,:,2,ci);
                        IMseriesMFIChNorm = reshape(IMseriesMFIChNorm,[],PeriodNumR(ci));
                        %                 IMseriesMFIChNorm = IMseriesMFIChNorm';
                        % ######## After correction ########
                        IMseriesCorMFI = reshape(IMseriesCor,[],PeriodNumR(ci));
                        %                 IMseriesCorMFI = IMseriesCorMFI';
                    end
                    
                    for ri = 2:size(ROIs,2)
                        TempValue = mean(IMseriesRspChS(mask_id{ri},:),1,'omitnan');
                        RspOfSession.ChS{i,ri-1} = cat(2,RspOfSession.ChS{i,ri-1},TempValue'); % Response, excited by 488 or 561 nm
                        TempValue = mean(IMseriesRspChNorm(mask_id{ri},:),1,'omitnan');
                        RspOfSession.ChNorm{i,ri-1} = cat(2,RspOfSession.ChNorm{i,ri-1},TempValue'); % Response, excited by 405 nm
                        TempValue = mean(IMseriesRspCor(mask_id{ri},:),'omitnan');
                        RspOfSession.Cor{i,ri-1} = cat(2,RspOfSession.Cor{i,ri-1},TempValue'); % Response, after correction
                        if MFIrecordFlag
                            TempValue = mean(IMseriesMFIChS(mask_id{ri},:),1,'omitnan');
                            MFIOfSession.ChS{i,ri-1} = cat(2,MFIOfSession.ChS{i,ri-1},TempValue'); % Average fluorescence intensities, excited by 488 or 561 nm
                            TempValue = mean(IMseriesMFIChNorm(mask_id{ri},:),1,'omitnan');
                            MFIOfSession.ChNorm{i,ri-1} = cat(2,MFIOfSession.ChNorm{i,ri-1},TempValue'); % Average fluorescence intensities, excited by 405 nm
                            TempValue = mean(IMseriesCorMFI(mask_id{ri},:),1,'omitnan');
                            MFIOfSession.Cor{i,ri-1} = cat(2,MFIOfSession.Cor{i,ri-1},TempValue'); % Average fluorescence intensities, after correction
                        end
                    end
                    if MFIrecordFlag
                        save([SaveFolder1,'\',SessionChTitle{fi,ci},'_RspOfSession.mat'],'RspOfSession','MFIOfSession','-v7.3');
                    else
                        save([SaveFolder1,'\',SessionChTitle{fi,ci},'_RspOfSession.mat'],'RspOfSession','-v7.3');
                    end
                    toc,
                    clear RspOfSession IMseriesRsp IMseriesCor IMseriesRspCor IMseriesRspChNorm IMseriesRspChS IMseriesMFIChS IMseriesMFIChNorm IMseriesCorMFI;
                end
                %% Analysis of response of each trial
                % response of each ROIs in each single trials
                % calculate response
                if StimRspFlag % Analyze images based on the stimulus information
                    tic,
                    StimTag = StimTagList{fi};   % Name of stimulus
                    RspOfSingleTrial.ChS = cell(StimModeNum,size(ROIs,2)-1); % Excited by 488 or 561 nm
                    RspOfSingleTrial.ChNorm = cell(StimModeNum,size(ROIs,2)-1); % Excited by 405 nm
                    RspOfSingleTrial.Cor = cell(StimModeNum,size(ROIs,2)-1); % After correction
                    Chfi = 1; % Excited by 488 or 561 nm
                    ChNormIdx = 2; % Excited by 405 nm
                    for i = 1:StimModeNum
                        if StimModeCount(i)>0
                            IMseriesTrialRspChS{i} = zeros([ImSz1,ImSz2,TrialFrames,StimModeCount(i)],'single');
                            IMseriesTrialRspChNorm(i) = IMseriesTrialRspChS(i);
                            IMseriesTrialRspCor(i) = IMseriesTrialRspChS(i);
                            IdxStim = find(ImStamp.Flag(:,1) == i);
                            IdxTrial = ImStamp.Note(IdxStim,2);
                            IdxStim = IdxStim-PeriodNumRange(1)+1;
                            IdxStart = IdxStartAll(IdxTrial);
                            IdxEnd = IdxEndAll(IdxTrial);
                            for ti = 1:StimModeCount(i) % Total number of trials for a specific type of stimulus
                                IMseriesTrialChS = single(IMseries(:,:,IdxStart(ti):IdxEnd(ti),Chfi,ci));
                                IMbasal = mean(IMseriesTrialChS(:,:,1:IdxStim(ti)-IdxStart(ti)-IdxCorLen,Chfi),3);
                                RspChS  = IMseriesTrialChS./IMbasal-1;
                                IMseriesTrialRspChS{i}(:,:,:,ti) = RspChS;
                                IMseriesTrialChNorm = single(IMseries(:,:,IdxStart(ti):IdxEnd(ti),ChNormIdx,ci));
                                IMbasal = mean(IMseriesTrialChNorm(:,:,1:IdxStim(ti)-IdxStart(ti)-IdxCorLen,Chfi),3);
                                RspChNorm = IMseriesTrialChNorm./IMbasal-1;
                                IMseriesTrialRspChNorm{i}(:,:,:,ti) = RspChNorm;
                                
                                switch CorrectMethodFlag
                                    case 1
                                        IMseriesTrialCor = IMseriesTrialChS./IMseriesTrialChNorm;
                                    case 2
                                        IMbasalfit = thetaIMk.*IMbasal+thetaIMb;
                                        IMseriesTrialCor = IMseriesTrialChS-thetaIMk.*IMseriesTrialChNorm-thetaIMb+IMbasalfit;
                                end
                                
                                IMbasal = mean(IMseriesTrialCor(:,:,1:IdxStim(ti)-IdxStart(ti)-IdxCorLen),3);
                                RspCor = IMseriesTrialCor./IMbasal-1;
                                IMseriesTrialRspCor{i}(:,:,:,ti) = RspCor;
                                % response of each ROIs in each single trials
                                RspChS = reshape(RspChS,[],TrialFrames);
                                RspChS = RspChS';
                                RspChNorm = reshape(RspChNorm,[],TrialFrames);
                                RspChNorm = RspChNorm';
                                RspCor = reshape(RspCor,[],TrialFrames);
                                RspCor = RspCor';
                                for ri = 2:size(ROIs,2)
                                    RspOfSingleTrial.ChS{i,ri-1} = cat(2,RspOfSingleTrial.ChS{i,ri-1},mean(RspChS(:,mask_id{ri}),2)); % 记录通道的反应
                                    RspOfSingleTrial.ChNorm{i,ri-1} = cat(2,RspOfSingleTrial.ChNorm{i,ri-1},mean(RspChNorm(:,mask_id{ri}),2)); % 校正通道反应
                                    RspOfSingleTrial.Cor{i,ri-1} = cat(2,RspOfSingleTrial.Cor{i,ri-1},mean(RspCor(:,mask_id{ri}),2)); % 比值的反应
                                end
                            end
                            for ri = 2:size(ROIs,2)
                                % Excited by 488 or 561 nm
                                Value = RspOfSingleTrial.ChS{i,ri-1}; % Response
                                ValueM = mean(Value,2); % Averaged response
                                ValueSEM = std(Value,0,2,'omitnan')/sqrt(size(Value,2)); % SEM
                                RspOfSingleTrial.ChS_M{i,ri-1} = cat(2,ValueM,ValueSEM);
                                % Excited by 405 nm
                                Value = RspOfSingleTrial.ChNorm{i,ri-1}; % Response
                                ValueM = mean(Value,2); % Averaged response
                                ValueSEM = std(Value,0,2,'omitnan')/sqrt(size(Value,2)); % SEM
                                RspOfSingleTrial.ChNorm_M{i,ri-1} = cat(2,ValueM,ValueSEM);
                                % After correction
                                Value = RspOfSingleTrial.Cor{i,ri-1}; % Response
                                ValueM = mean(Value,2); % Averaged response
                                ValueSEM = std(Value,0,2,'omitnan')/sqrt(size(Value,2)); % SEM
                                RspOfSingleTrial.Cor_M{i,ri-1} = cat(2,ValueM,ValueSEM);
                            end
                            
                            % averaged response of each pixels
                            %                     if SaveSingleChFlag(ci)
                            IMseriesTrialRspChS_M = mean(IMseriesTrialRspChS{i},4);
                            FileName = fullfile(SaveFolder1,[SessionChTitle{fi,ci},'_',ChTag{ci,Chfi},'_',StimTag{i},'_TrialAveRsp.tif']);
                            F_WriteBigTiff_Ver1(FileName,pixelsize,compression,IMseriesTrialRspChS_M);
                            IMseriesTrialRspChNorm_M = mean(IMseriesTrialRspChNorm{i},4);
                            FileName = fullfile(SaveFolder1,[SessionChTitle{fi,ci},'_',ChTag{ci,ChNormIdx},'_',StimTag{i},'_TrialAveRsp.tif']);
                            F_WriteBigTiff_Ver1(FileName,pixelsize,compression,IMseriesTrialRspChNorm_M);
                            %                     else
                            %                         disp('No saving response data as tiff file based on singel channel (Trial average)');
                            %                     end
                            IMseriesTrialRspCor_M = mean(IMseriesTrialRspCor{i},4);
                            FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},ChTag{ci,ChNormIdx},'cor_',StimTag{i},'_TrialAveRsp','.tif']);% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
                            F_WriteBigTiff_Ver1(FileName,pixelsize,compression,IMseriesTrialRspCor_M);
                        end
                    end
                    RspOfSingleTrial.TimeInSec = TrialTimeFrameRcd;
                    RspOfSingleTrial.TimeInMin = RspOfSingleTrial.TimeInSec/60;
                    RspOfSingleTrial.TimeInHrs = RspOfSingleTrial.TimeInSec/3600;
                    FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},ChTag{ci,ChNormIdx},'cor_TrialRsp.mat']);
                    save(FileName,'IMseriesTrialRspCor','-v7.3');
                    save([SaveFolder1,'\',SessionChTitle{fi,ci},'_IMseriesOfSingleTrial.mat'],'RspOfSingleTrial','-v7.3');
                    clear IMseriesTrialRspCor RspOfSingleTrial IMseriesTrialRspChNorm IMseriesTrialRspChS;
                end
            end
            if exist('BackgroundRcd','var')
                save([SaveFolder1,'\',SessionChTitle{fi,1}(1:end-1),'BackgroundRcd.mat'],'BackgroundRcd','-v7.3');
            end
            if exist('TrialTimeFrameRcd','var')
                save([SaveFolder1,'\',SessionChTitle{fi,1}(1:end-1),'Time Of Each Frame.mat'],'TrialTimeFrameRcd','-v7.3');
                clear TrialTimeFrameRcd;
            end
            if exist('shiftsRcd','var')
                save([SaveFolder1,'\',SessionChTitle{fi,1}(1:end-1),'ShiftsRcd Of Each Frame.mat'],'shiftsRcd','-v7.3');
                clear shiftsRcd;
            end
            clear IMseries shiftsRcd targetsName; % 释放内存
    end
    tEnd = toc(tStart0);
    disp('Finished!'),
    disp(['Total time cost: ',num2str(tEnd/60),' min for ',num2str(SessionNum),' session;',newline,'For each session: ',num2str(tEnd/60/SessionNum),' min']),
    cd(MfileDir);
    %% =======================================================================================
    %% Function scanDir for scanning files or folders with specific string
    % scanDir, a function for scanning files or folders that match a specific name
    function TargetList = scanDir(CurrentPath,TargetLabel,TargetType,MatchPattern)
    % Serach folders and files under CurrentPath
    % CurrentPath = mainPath;
    files = dir(CurrentPath);
    len = length(files);
    TargetList = {};
    index = 1;
    for i = 1:len
        % Skip the . and .. directories
        if (strcmp(files(i).name, '.') == 1) ...
                || (strcmp(files(i).name, '..') == 1)
            continue;
        end
        % Recursively call the function to traverse the folders in the current directory (too deep recursion may lead to errors)
        if files(i).isdir == 1
            CurrentPath1 = [CurrentPath, '\', files(i).name];
            disp(CurrentPath1);
            tmpName = scanDir(CurrentPath1,TargetLabel,TargetType,MatchPattern);
            for j = 1 : length(tmpName)
                TargetList{index} = tmpName{j};
                index = index + 1;
            end
        end
        switch TargetType
            case 1
                if MatchPattern
                    if files(i).isdir == 1 && strcmp(files(i).name,TargetLabel)
                        TargetList{index} = fullfile(CurrentPath, '\', files(i).name);
                        index = index + 1;
                    end
                else
                    if files(i).isdir == 1 && ~isempty(strfind(files(i).name,TargetLabel))
                        TargetList{index} = fullfile(CurrentPath, '\', files(i).name);
                        index = index + 1;
                    end
                end
            case 2
                if MatchPattern
                    if files(i).isdir == 0 && strcmp(files(i).name,TargetLabel)
                        TargetList{index} = fullfile(CurrentPath, '\', files(i).name);
                        index = index + 1;
                    end
                else
                    if files(i).isdir == 0 && ~isempty(strfind(files(i).name,TargetLabel))
                        TargetList{index} = fullfile(CurrentPath, '\', files(i).name);
                        index = index + 1;
                    end
                end
            case 3
                if MatchPattern
                    if strcmp(files(i).name,TargetLabel)
                        TargetList{index} = fullfile(CurrentPath, '\', files(i).name);
                        index = index + 1;
                    end
                else
                    if ~isempty(strfind(files(i).name,TargetLabel))
                        TargetList{index} = fullfile(CurrentPath, '\', files(i).name);
                        index = index + 1;
                    end
                end
        end
    end
    end
    
    %% =======================================================================================
    %% F_WriteBigTiff_Ver1
    % For outputing large TIFF files, capable of handling files larger than 4GB
    % Fei Deng,20210808
    function F_WriteBigTiff_Ver1(FileName,pixelsize,compression,ImStack)
    obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
    [ImSz1 ImSz2 ImSz3] = size(ImStack);
    for i = 1:ImSz3
        img = ImStack(:,:,i);
        img = img'; % Transpose the image before saving it (WriteIMG leads to inconsistent after saving)
        WriteIMG(obj,img);
    end
    close(obj);
    end
    %% =======================================================================================
    %% F_StimModeImport_Ver0,used for read/import StimMode record from the txt file
    % A function to read the "StimMode" records from a text file. Fei Deng,20210620
    function StimMode = F_StimModeImport_Ver0(fullFileName,SplitTag)
    % SplitTag = 'StimMode: ';
    fid = fopen(fullFileName,'rt');
    % count the number of lines in the text file
    txtRow = 0;
    while ~feof(fid)
        txtRow = txtRow+sum(fread(fid,10000,'*char')==char(10));
    end
    fclose(fid);
    disp(['There are ',num2str(txtRow),' lines in txt file, importing data:']);
    % read the raw data from the text file
    recordV = cell(txtRow,2); % Preallocate storage space to improve computation speed
    fidin = fopen(fullFileName);
    i = 0;
    % SplitTag = '2';
    while ~feof(fidin)
        tline = fgetl(fidin); %读取每行文本
        
        %     if isempty(tline)
        %         continue;
        %     end
        %     if tline(1)~= SplitTag
        %         continue;
        %     end
        i = i+1;
        if mod(i,10000) == 0
            disp(i);
        end
        %  m = split(tline,'->'); % Separate the columns of each row
        m = split(tline,SplitTag); % Separate the columns of each row
        recordV(i,1:size(m,1)) = m';
    end
    fclose(fidin);
    StimModeIdx = find(~cellfun(@isempty,recordV(:,2)));
    for i = 1:size(StimModeIdx,1)
        StrTemp = recordV{StimModeIdx(i),2};
        StimMode(i) = str2double(StrTemp(1:find(StrTemp==',')-1));
    end
    StimMode = StimMode';
    disp('Importing StimMode data finished');
    end
    %% =======================================================================================
    %% F_RefDarkImport_Ver0
    %% ## Fei Deng,20220117,used for import Dark images then generate offset images
    function RefDark = F_RefDarkImport_Ver0(TargetPathD,filetype,pixR,pixC,CameraNum,pixelSzBin1,ImBinning,compression)
    disp('Looking for dark images in folder:');
    disp(TargetPathD);
    RefDark = uint16([]);
    try
        DarkDir = fullfile(TargetPathD,'Dark');
        cd(DarkDir);
        DarkList = dir('*Dark*');
        if size(DarkList,1) == CameraNum
            for ci = 1:CameraNum
                targetsName = DarkList(ci).name;
                if strfind(targetsName,filetype)
                    img = imread(targetsName);
                    RefDark(:,:,ci) = img;
                    figure,imshow(img,[]);
                    title(DarkList(ci).name,'Interpreter','none');
                    disp([targetsName,' exist']);
                else
                    FolderTemp = fullfile(DarkList(ci).folder,targetsName,'Default');
                    cd(FolderTemp);
                    targetPack = dir(['*.',filetype]);
                    targetsName = {targetPack(:).name};
                    targetsName = targetsName';
                    IMseries = zeros(round([pixR,pixC,size(targetsName,1)]),'uint16');
                    for Prdi = 1:size(IMseries,3)
                        FileName = targetsName{Prdi};
                        IMseries(:,:,Prdi) = imread(FileName);
                    end
                    FileName = fullfile(DarkDir,[DarkList(ci).name,'.tif']);
                    img = uint16(mean(IMseries,3));
                    RefDark(:,:,ci) = img;
                    figure,imshow(img,[]);
                    title(DarkList(ci).name,'Interpreter','none');
                    obj = Fast_BigTiff_Write(FileName,pixelSzBin1/ImBinning,compression);
                    WriteIMG(obj,img'); % save averaged image
                    close(obj);
                    cd(DarkDir);
                    rmdir(DarkList(ci).name, 's');
                end
            end
            disp([FileName,' generated and raw folder deleted.']);
        else
            f = errordlg('Dark image inconsistent with camera number','Dark image Error');
        end
    catch
        disp('No Dark images!');
    end
    end
    
    %% =======================================================================================
    %% F_BasalIM_Ver1
    %% ## Fei Deng,20220328,used for generating IMseries3DBasalIM image
    % Fei Deng,20220413,automatic detection of the input data type
    %% =======================================================================================
    function [BasalIM,BasalIM_Thd] = F_BasalIM_Ver1(IMseries3D,BasalFlagOfSession,LowThd,IdxBaseFrame)
    className = class(IMseries3D); % automatic detection of data type
    [ImSz1,ImSz2,ImSz3] = size(IMseries3D);
    BasalIM = zeros(ImSz1,ImSz2);
    BasalIM_Thd = BasalIM;
    % IMseries3DBasalIM = zeros(ImSz1,ImSz2),ImSz3,'uint16');
    switch BasalFlagOfSession
        case 3
            tic,
            for coli = 1:ImSz2
                %                                     coli = 150
                img = squeeze(IMseries3D(:,coli,:));
                %                             figure,imshow(img,[]);
                imgSort = sort(img,2);
                %                             figure,imshow(imgSort,[]);
                BasalIM(:,coli) = mean(imgSort(:,1:round(LowThd*ImSz3)),2);
                BasalIM_Thd(:,coli) = imgSort(:,round(LowThd*ImSz3));
            end
        otherwise
            BasalIM = squeeze(mean(IMseries3D(:,:,IdxBaseFrame),3));
    end
    BasalIM = eval([className,'(BasalIM)']); % Change the data type to match the input data type
    BasalIM_Thd = eval([className,'(BasalIM_Thd)']); % Change the data type to match the input data type
    % figure,imshow(BasalIM,[]);
    % title('Basal image');
    % close;
    end
