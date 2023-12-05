%% MesoscopicAnalysis_rawdata
% Running time:20230925
%   ======================================================================================
% MesoscopicAnalysis_rawdata,20230925，Fei Deng @ Yulong Li lab edited
%   =======================================================================================
%%
close all,
clear,clc;

%% Parameter setting
MotionCorFlag = 1; % whether do Motion correction, 1-Yes, 0-No 
RegisterFlag = 1; % Whether register images from different cameras, 1-Yes, 0-No 
BaselineCorFlag = 0; % Whether correct the baseline，1-use functions in Matlab toolbox，e.g. exp2; 2-use beads function in Matlab toolbox；3-use smooth function (very slow)
SpectraUnmixFlag = 1; % Whether do spectra unmixing, 1-Yes, 0-No 
CorrectMethodFlag = 2; % Methods for correction based on 405 nm light excited images, 1-ratio; 2-fit then substract
WholeSessionRspFlag = 1; % Whether analyze response if the whole session, 1-Yes, 0-No 
StimRspFlag = 0; % Whether analyze response according different stimuli
BasalFlagOfSession = 1;  % Methods for generating the baseline image for the whole session: 1-Assignment frame；2-Based on the stimulus marker；3. X% low value for each pixel the whole session; 4-The union set of baseline from the begin of session and each trial; 5-based on other labelled baseline, such as REM during sleep wake cycles
IdxBaseList = ones(20,2); % 每个session实际输入并且处理的图片帧数范围，初始置最小值避免报错，第一列为开始帧数，第二列为结束帧数
IdxBaseList(:,2) = 9999999999; % 初始置最大值避免报错
IdxBaseList(1,:) = floor([1,35]); % 根据实际情况手动设定指定session的baseline对应帧数,按照每个通道的图片编号
forceSaveFlag = 1; % 是否以二进制文件保存IMseries
SaveInSubFolderFlag = 0; % 是否将不同session的结果分别存放在不同文件夹
MFIrecordFlag = 1; % 是否分析整个session的平均荧光强度,0-不分析；
SaveSingleChFlag = [1,1]; % 是否保存每个相机（绿色、红色）的单个通道的反应数据（未校正的反应数据）,0-不保存，1-保存；
% mainPath = 'K:\Mouse in vivo\20210104_L21\SleepWake'; % 输入查找文件夹范围
TargetLabel = 'Default'; % 输入目标文件夹/文件特征字符串
MatchPattern = 1; % scanDir函数是否需要字符全匹配
TargetType = 1; % 1--文件夹，2--文件，3--文件夹和文件
filetype = 'tif';% 图片格式
% #################### 成像信息 ####################
CameraNum = 2; % 使用的相机数目
ChStatus = [1,1,1]; % 488nm, 561nm, 405nm for imaging
ChSim = [1,1,0]; % 488nm, 561nm, 405nm 同时成像的通道
ChNorm = 3; % Normalized channel, 1-488nm, 2-561nm, 3-405nm
CameraTag = {'G','R','F'};
ChName = {'Ex488','Ex561','Ex405'};
ChNameTag = {'B','G','V'};  % 激发光缩写
ChColor = {'g','r','b'};
fixCh = 1; % 用于配准的固定通道编号
ExposurePrecision = -2; % 根据实际曝光时间精度调整,如0.035s对应精度为-3，如0.04s对应精度为-2
ImSampleModeNum = 1; % 成像采样频率种数
Cgr = 0.119602778; % 绿色荧光信号串至红色通道比例，基于5-HT3.0实测值0.119602778
% Cgr = 0.098192171; % 绿色荧光信号串至红色通道比例，基于eCB2.0实测校正值 Cgr = 0.098192171;
% Cgr = 0; % 绿色荧光信号串至红色通道比例，基于eCB2.0实测校正值 Cgr = 0.098192171; 此处由于成像设置，理论上绿色sensor不会串到红色通道
% Cgr = 0.073863753; % 绿色荧光信号串至红色通道比例，基于EGFP实测校正值 Cgr = 0.073863753
Crg = 0.001328723; % 红色荧光信号串至绿色通道比例，基于jRGECO1a = 0.001328723,
% Crg = 2.69021E-07; % 红色荧光信号串至绿色通道比例，基于r5-HT2.0 = 2.69021E-07
ChName = ChName(find(ChStatus==1));
ChColor = ChColor(find(ChStatus==1));
fixCh = sum(ChStatus(1:fixCh));
ChNumAll = length(ChName); % 激发光总数
ChSimNum = sum(ChSim);
ChNum = ChNumAll-(ChSimNum-1); % 记录的每周期激发光激发次数
ChNorm0 = min([ChNum,ChNorm]); % 20230619修改此处为ChNorm0，以便和后续的ChNorm = 2区分开，避免后续赋值后被改变而报错
switch ChNum
    case 2
        ChTag0 = {'B','V';'G','V';};
    case 3
        ChTag0 = repmat(ChNameTag,CameraNum,1);
end

IMscale = 0.7; % Resize图片系数  IMscale = 0.7;
PeriodNumRangeList = ones(10,2); % 每个session实际输入并且处理的图片帧数范围，初始置最小值避免报错
PeriodNumRangeList(:,2) = 9999999999; % 初始置最大值避免报错
% PeriodNumRangeList(1,2) =200; % 根据实际需要手动修改PeriodNumRangeList的某些值
ImBinning = 4; % 实际成像时的binning数值  ImBinning = 3;
pixelSzBin1 = 2211.358; % 成像时binning为1X1时的标尺，pixels / cm， 20211031测定校准标尺
pixelsize0 = pixelSzBin1/ImBinning;
pixelsize = pixelSzBin1/ImBinning*IMscale;
compression = 0;
% BgToAutoF = [2,2,2;2,2,2]; % 不表达任何探针的小鼠的皮层自发荧光和血管的荧光比值，每行为不同相机通道，每列为不同激发光通道

% #################### 刺激信息（如无特定刺激，可以不用修改） ####################
TrialFastTime = 120;  % 每个trial的高速采样时长（s）
StimT = 10;  % 最大刺激时长（s）
StimMaxT = StimT*1.2;  % 最大刺激时长（s）
StimTagList = {{'Base_20Hz 1s','Base_20Hz 10s','Base_50Hz 10s','GBR_20Hz 1s','GBR_20Hz 10s','GBR_50Hz 10s'}};% 刺激名称，每个session一行
BaseDurList = ones(18,1)*10; %每个session的每个trial baseline计算时使用时间(s)
IdxCorLen = 2; % 为避免光激活时基线校正使用的数据受光激活激发光的影响，所选数据点尽量前移，前移IdxCorLen个数据点

% ResultTag = {'Miss','Hit','CR','FA'};

% #################### 固定参数（通常无需修改） ####################
BlockSz = 10000; % 分组处理数据，每组帧数
MfileDir = pwd;
addpath(genpath(MfileDir));
ParentFolder = fileparts(MfileDir);
mfileNameR = mfilename; % 当前m文件文件名
rawdataName = mfileNameR(strfind(mfileNameR,'rawdata'):end);
mainPath = fullfile(ParentFolder,rawdataName);
ImCorRatioFile = fullfile(MfileDir,'System Setting\CorRatio_20211017.mat'); % 20211017至20220613之间的数据使用此校正文件
% ImCorRatioFile = fullfile(MfileDir,'CorRatio_20220613.mat'); % 20220613之后的数据使用此校正文件
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

%% 共用图片信息读取
cd(TargetList{1});
targetPack = dir(['*.',filetype]);
IM1st = imread(targetPack(1).name);
[pixR,pixC] = size(IM1st); % 原始图片像素尺寸
% 读取offset图片
TargetPath = mainPath;
RefDarkA = uint16(zeros(pixR,pixC,CameraNum));
RefDarkA = F_RefDarkImport_Ver0(TargetPath,filetype,pixR,pixC,CameraNum,pixelSzBin1,ImBinning,compression);
%% 处理不同session
for fi = 1:SessionNum
    close all,
    TargetPath = TargetList{1+(fi-1)*CameraNum}; % 为方便记录，所有原始数据按照session放置
    %     TargetPath = TargetList{fi};
    [filepath,SessionTitle{fi,1},ext] = fileparts(fileparts(TargetPath));
    if ~isempty(ext) % 避免文件夹名称含小数点导致的fileparts函数识别出错的情况
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
    % 读取offset图片
    TargetPathD = fileparts(fileparts(TargetPath));
    RefDark = F_RefDarkImport_Ver0(TargetPathD,filetype,pixR,pixC,CameraNum,pixelSzBin1,ImBinning,compression);
    if isempty(RefDark)
        if isempty(RefDarkA)
            disp(['No Dark images for ', SessionTitle{fi,1},'!']);
            RefDark = uint16(zeros(pixR,pixC,CameraNum));
        else
            RefDark = RefDarkA; % 使用公用dark images
        end
    end
    %% Read Spike2 recording
    StimRcdPath = fullfile(filepath,'StimRcd'); % Path of StimRcd folder
    cd(StimRcdPath);
    Spk2RecordFile = dir('*.mat');
    Spk2RecordFile = Spk2RecordFile.name; % Spike2 recording task performance
    load(Spk2RecordFile);
    % ##########################################################################
    % 针对spike2记录的imaging数据点多出一个的debug，删除最后多出来的这个点
    if  sfi == 999
        Imaging.level = Imaging.level(1:end-1);
        Imaging.times = Imaging.times(1:end-1);
    end
    % #####################################################################
    ImagingInfo = cat(2,Imaging.times,double(Imaging.level));
    % 计算成像采样频率
    ImageStateLast = [ImagingInfo(2:end,1)-ImagingInfo(1:end-1,1)]; % 每一状态持续时间,从该状态（0或者1）改变开始计时
    ImageStateLastR = roundn(ImageStateLast,ExposurePrecision); %根据实际曝光时间精度调整
    %     ImageStateLastR = roundn(ImageStateLast,-4);
    ImageStateLastRFreq = tabulate(ImageStateLastR);
    ImageStateLastRFreq = sortrows(ImageStateLastRFreq,3,'descend');
    % 针对spike2记录的帧数比实际成像多帧的debug，找出多出来的数据点为BugIdx1至BugIdx2
    BugDurIdx = find(ImageStateLastRFreq(:,3)<0.01); % 出现频率很低的时间间隔点应该异常
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
            BugIdx = BugIdx+1; % 对应ImagingInfo的行数
        else
            BugFlag = 0;
        end
        if BugFlag
            figure,plot(ImageStateLastR(BugIdx(1)-1:max(BugIdx(end))+1));
            ImagingInfo0 = ImagingInfo;
            ImageStateLastR0 = ImageStateLastR;
            %         ImagingInfo = ImagingInfo0;
            %         figure,plot(ImagingInfo0(3000:3022,1),ImagingInfo0(3000:3022,2));
            ImagingInfo(BugIdx,:) = []; % 针对spike2记录的帧数比实际成像多帧debug，多出来的数据点被删除
            %             ImagingInfo(BugIdx(1:end-1),:) = []; % 针对spike2记录的帧数比实际成像多帧debug，多出来的数据点被删除
            %          ImagingInfo(BugIdx(1:end)+1,:) = []; % 针对spike2记录的帧数比实际成像多1帧debug，多出来的数据点为66263，66264
            ImageStateLast = [ImagingInfo(2:end,1)-ImagingInfo(1:end-1,1)]; % 每一状态持续时间,从该状态（0或者1）改变开始计时
            ImageStateLastR = roundn(ImageStateLast,-2); %根据实际曝光时间精度调整
            ImageStateLast = [ImagingInfo(2:end,1)-ImagingInfo(1:end-1,1)]; % 每一状态持续时间,从该状态（0或者1）改变开始计时
            ImageStateLastR = roundn(ImageStateLast,-2); %根据实际曝光时间精度调整
            %     ImageStateLastR = roundn(ImageStateLast,-4);
            ImageStateLastRFreq = tabulate(ImageStateLastR);
            ImageStateLastRFreq = sortrows(ImageStateLastRFreq,3,'descend');
        end
        BugDurIdx = find(ImageStateLastRFreq(:,3)<0.01); % 出现频率很低的时间间隔点应该异常
    end
    % 成像信息
    ExTotalNum = size(ImagingInfo,1)/2; % 总共给出的激发光脉冲数目
    ImagingTimeM = mean(reshape(ImagingInfo(:,1),ChNum*2,[])); % 根据通道数计算成像时刻
    BlockFrameNum = size(ImagingTimeM,2);
    %     ImSampleModeNum = size(ImageStateLastRFreq,1)-2;
    Period_T = zeros(1,ImSampleModeNum);
    if ChNum==1
        for i =  1:ImSampleModeNum
            Period_T(i) = ImageStateLastRFreq(1,1)+ImageStateLastRFreq(2+i,1); % 每次成像周期（含成像和间隔时间）,单位s
        end
    elseif ChNum>1
        for i =  1:ImSampleModeNum
            Period_T(i) = ImageStateLastRFreq(1,1)*ChNum+ImageStateLastRFreq(2,1)*(ChNum-1)+ImageStateLastRFreq(2+i,1); % 每次成像周期（含成像和间隔时间）,单位s
        end
    end
    Period_T = sort(Period_T,'ascend');
    Image_f = 1./Period_T; % 成像频率/Hz
    
    if StimRspFlag % 读取刺激信息并分析整理
        StimInfo0 = cat(2,Stim.times,double(Stim.level));
%         Stim.times(Stim.length-1:Stim.length) = []; % 针对结尾处多出了一个高电平记录进行debug
%         Stim.level(Stim.length-1:Stim.length) = []; % 针对结尾处多出了一个高电平记录进行debug
        TagColor = {'black','red'}; % red-StimInfo, black-ImagingInfo
        % 计算刺激数目并找到刺激起点,生成仅记录每个刺激起止时刻的StimInfo
        StimStateLast = [Stim.times(2:end)-Stim.times(1:end-1)]; % 每一状态持续时间,从该状态（0或者1）改变开始计时
        StimEndIdx = [find(StimStateLast>StimMaxT);size(Stim.times,1)]; % 每个trial的刺激结束时刻，在原始刺激电平记录结果中的索引
        TrialNum = size(StimEndIdx,1);
        StimStartIdx = [1;StimEndIdx(1:end-1)+1]; % 每个trial的刺激开始时刻，在原始刺激电平记录结果中的索引
        StimInfo = zeros(TrialNum*2,2);
        StimInfo(1:2:TrialNum*2-1,1) = Stim.times(StimStartIdx); % 每个trial的刺激开始时刻
        StimInfo(1:2:TrialNum*2-1,2) = 1;
        %     StimInfo(2:2:TrialNum*2,1) = Stim.times(StimEndIdx); % 每个trial的刺激结束时刻
        MultiPulseIdx = find(StimEndIdx-StimStartIdx>1); % 查找含有多个刺激信号脉冲的trial编号
        if MultiPulseIdx % 判断每个trial是否有多个刺激信号脉冲
            %              if size(StimInfo0,1)>size(StimInfo,1) % 判断每个trial是否有多个刺激信号脉冲
            StimInfo(MultiPulseIdx'*2,1) = Stim.times(StimEndIdx(MultiPulseIdx))+(Stim.times(StimStartIdx(MultiPulseIdx)+2)-Stim.times(StimStartIdx(MultiPulseIdx)+1)); % 每个trial的刺激结束时刻，加上非高电平脉冲的时间以作校正
        else
            StimInfo(2:2:TrialNum*2,1) = Stim.times(StimEndIdx); % 每个trial的刺激结束时刻，每个trial只有一个脉冲时，不作校正
        end
        %         try
        %             StimInfo(2:2:TrialNum*2,1) = Stim.times(StimEndIdx)+(Stim.times(StimStartIdx+2)-Stim.times(StimStartIdx+1)); % 每个trial的刺激结束时刻，加上非高电平脉冲的时间以作校正
        %         catch exception
        %             StimInfo(2:2:TrialNum*2,1) = Stim.times(StimEndIdx); % 每个trial的刺激结束时刻，每个trial只有一个脉冲时，不作校正
        %         end
        StimInfo(2:2:TrialNum*2,2) = 0;
        
        %     if StimInfo(2,1)-StimInfo(1,1)>1.5*StimDur % 第一个矩形脉冲时间偏离刺激持续时间过大时，应该为arduino启动时的标志，舍去
        %         StimInfo = StimInfo(3:end,:);
        %     end
        
        %     TrialDur = 360;
        if TrialNum==1
            TrialDur = 360; % 时间可能需根据实际情况修改？
        else
            TrialDur = (StimInfo(end-1,1)-StimInfo(1,1))/(TrialNum-1); %每个trial平均总时间(s)
        end
        StimDur = StimInfo(2:2:end,1)-StimInfo(1:2:end-1,1); %每个trial设定的刺激时长，单位s
        TrialDur = round(TrialDur); %设定的理论值为整数
        StimDur = round(StimDur); %设定的理论值为整数
        BaseDur = BaseDurList(fi,1); %分析时设定的每个session的每个trial baseline时间(s)
        RestDur = TrialDur-BaseDur; %每个刺激起始至结束后时间(s)
        
        %     StimInfo = StimInfo(1:TrialNum*2,:); % 删除arduino未记录的trial
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
        % 标注成像时刻
        yTemp = ones(1,BlockFrameNum)*0.85;
        text(ImagingTimeM,yTemp,num2cell(1:BlockFrameNum)); % 标记imaging period编号
        StimTimeM = mean(reshape(StimInfo(:,1),2,[]));
        yTemp = ones(1,TrialNum)*0.95;
        text(StimTimeM,yTemp,num2cell(1:TrialNum),'Color','red'); % 标记imaging period编号
        legend({'Imaging','Stim'});
        xlabel('Time (s)');
        saveas(gcf,[SaveFolder1,'\',SessionTitle{fi,1},'_Stim&Imaging Info.fig']);
        saveas(gcf,[SaveFolder1,'\',SessionTitle{fi,1},'_Stim&Imaging Info.jpg']);
        %% State labelling on images
        % 结构体ImStamp记录每一帧图像对应的不同时间的标签及事件内部编号,按照列排序分别为Lick,Sound,Window,Actor,Result
        %     ImStamp.Flag = int8(zeros(BlockFrameNum,size(TaskData,2)));
        %     ImStamp.Note = int8(zeros(BlockFrameNum,size(TaskData,2)));
        ImStamp.Flag = zeros(BlockFrameNum,1)*nan;
        ImStamp.Note = ImStamp.Flag;
        % ####### 图片贴标签方案1：将每种事件标签贴至该事件发生后的最近一帧图像 #######
        % Stimuli标签
        ArdRcdFile = [Spk2RecordFile(1:end-4), '_ArdRcd.txt'];
        fullFileName = fullfile(StimRcdPath,ArdRcdFile);
        if exist(ArdRcdFile)
            SplitTag = 'StimMode: ';
            %             SplitTag = 'StimMode '; % 针对记录文件格式异常修改
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
            CurImIdx = find(ImagingTimeM>StimInfo(2*i-1,1),1); % 找到大于每个刺激起始时间的图片编号
            if CurImIdx
                ImStamp.Flag(CurImIdx,1) = StimMode(i,1); % 该类型事件种类编号,即刺激模式编号
                StimModeCount(StimMode(i,1)) = StimModeCount(StimMode(i,1))+1;
                ImStamp.Note(CurImIdx,1) = StimModeCount(StimMode(i,1)); % 该类型事件序号，即属于该事件/刺激的第几个trial
                ImStamp.Note(CurImIdx,2) = i; % 总trial数编号
            end
        end
        IdxStimAll = find(ImStamp.Note(:,2)>0);
        %     TrialFrames = round(mean(diff(IdxStimAll)));  % 每个trial的平均采样帧数
        TrialFrames = min(diff(IdxStimAll));  % 每个trial的最小采样帧数
        %         TrialFrames = max(diff(IdxStimAll));  % 每个trial的最大采样帧数
        %         TrialFrames = round(mean(diff(IdxStimAll)));  % 每个trial的平均采样帧数
        IdxStartAll = IdxStimAll-BaseDur*Image_f(1);
        IdxEndAll = IdxStartAll+TrialFrames-1;
        if length(Period_T)>1
            TrialFastFrames = TrialFastTime*Image_f(1); % 每个trial的平均快速采样帧数
            IdxFastSampleEndAll = IdxStartAll+TrialFastFrames-1;
            TrialTimeFrameRcd = cat(2,(1:TrialFastFrames)*Period_T(1),TrialFastFrames*Period_T(1)+(1:TrialFrames-TrialFastFrames)*Period_T(2)); % 每一个trial的每帧图像的记录相对时间
        else
            TrialTimeFrameRcd = (1:TrialFrames)*Period_T;
        end
        TrialTimeFrameRcd = TrialTimeFrameRcd-BaseDur;  % 以刺激时刻为0时刻，刺激后的第一帧为第一次大于0的时刻
        TrialTimeFrameRcd = TrialTimeFrameRcd';
        save([SaveFolder1,'\',SessionTitle{fi,1},'_Stim&Imaging Info.mat'],'-v7.3');
    end
    
    %% Import imaging info
    ChAnalyse = zeros(CameraNum,ChNum); % 488nm, 561nm, 405nm for analyse
    for ci = 1:CameraNum % 同一个session根据每个相机记录的数据分别处理，CameraNum为相机数目
        switch ChSimNum
            case 1
                ChAnalyse(ci,[ci,ChNorm0]) = 1; % 488nm, 561nm, 405nm for analyse, for camera ci
            case 2
                ChAnalyse(ci,[1,ChNorm0]) = 1; % 488nm, 561nm, 405nm for analyse, for camera ci
        end
        ChAnlsNum(ci) = mean(sum(ChAnalyse(ci,:),2)); % channel number for analyse
        TargetPath = TargetList{ci+(fi-1)*CameraNum}; % 为方便记录，每个相机记录文件夹名称格式为XXX_G_1，所有原始数据按照session放置
        %           TargetPath = TargetList{fi+(ci-1)*SessionNum};
        [filepath,SessionChTitle{fi,ci},ext] = fileparts(fileparts(TargetPath));
        if ~isempty(ext) % 避免文件夹名称含小数点导致的fileparts函数识别出错的情况
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
        % #########读取原始图片#########
        filesNO(ci) = size(targetPack,1);
        PeriodNum(ci) = floor(filesNO(ci)/ChNum); %每个相机的成像周期数
        % 删除无用图片（针对红绿通道的561和488nm激发光不同时激发的情况
        if filesNO(ci) == ExTotalNum && ~isempty(find(ChAnalyse(ci,:)==0))
            %         if filesNO(ci) >70000 %手动修改判断条件
            tic,
            disp(['Deleting invalid images......']);
            for Chfi = find(ChAnalyse(ci,:)==0) %1:ChNum  %读取同一相机记录文件夹下的图片
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
            % 重新获取删除无用图片后的图片文件列表
            targetPack = dir(['*.',filetype]);
            targetsName(:,ci) = ({targetPack(:).name})';
            filesNO(ci) = size(targetPack,1);
        else
            disp('No images for deleting');
            targetsName(:,ci) = targetsNameTemp;
        end
        PeriodNum(ci) = floor(filesNO(ci)/ChAnlsNum(ci)); %每个相机的成像周期数
        PeriodNumRange = PeriodNumRangeList(fi,1):min([PeriodNumRangeList(fi,2),PeriodNum(ci)]); % 实际处理的图片帧数范围
        PeriodNumR(ci) = length(PeriodNumRange);
        % 针对Spike2 记录时长短于实际成像时长的bug增加该段用于校正
        if BlockFrameNum<PeriodNum(ci)
            disp('Bug! Spike2 记录时长短于实际成像时长');
            ImagingTimeM0 = ImagingTimeM;
            ImagingTimeM(BlockFrameNum+1:PeriodNum(ci)) = ImagingTimeM(end)+max(Period_T)*(1:PeriodNum(ci)-BlockFrameNum);
            disp('成像时刻已校正！');
            BlockFrameNum = size(ImagingTimeM,2);
        end
    end
    ImagingTimeM0 = ImagingTimeM;
    ImagingTimeM = ImagingTimeM(PeriodNumRange);
    
    if StimRspFlag % 根据待处理图像帧数范围调整IdxStartAll和IdxEndAll的值
        IdxStartAll0 = IdxStartAll;
        IdxEndAll0 = IdxEndAll;
        IdxStartAll = IdxStartAll-PeriodNumRange(1)+1;
        IdxEndAll = IdxEndAll-PeriodNumRange(1)+1;
        BadIdx = find(IdxEndAll>PeriodNumRange(end));
        IdxStartAll(BadIdx) = [];
        IdxEndAll(BadIdx) = [];
    end
    %% 对不同相机记录结果进行处理
    BlockNum = ceil(PeriodNumR(1)/BlockSz); % 分组读取数据的组数
    IMseriesTave =  zeros(pixR,pixC,CameraNum,ChAnlsNum(1),BlockNum,'uint16');
    IMseriesTmin = IMseriesTave;
    IMseriesTmax = IMseriesTave;
    IMseriesRawAve = zeros(pixR,pixC,CameraNum,ChAnlsNum(1),'uint16');
    IMseriesRawMin = IMseriesRawAve;
    IMseriesAve  = zeros([ceil([pixR*IMscale,pixC*IMscale]),CameraNum,ChAnlsNum(1)],'uint16');
    
    % 生成 IMseries
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
        IMseries = F_load_IMseries_from_datbin_file_Ver1(filename);  % 从已暂存的二进制文件中快速读取
    else
        baseCorFinishFlag = 0;
        IMseries = zeros([ceil([pixR*IMscale,pixC*IMscale]),PeriodNumR(ci),ChAnlsNum(ci),CameraNum],'uint16');
        shiftsRcd = [];
        tStart = tic;
        for ci = 1:CameraNum
            TempTitle = SessionChTitle{fi,ci};
            %         Idx = find(TempTitle=='_',1,'last');
            TargetPath = TargetList{ci+(fi-1)*CameraNum}; % 为方便记录，每个相机记录文件夹名称格式为XXX_G_1，所有原始数据按照session放置
            cd(TargetPath);
            disp(['Importing data folder: ',TargetPath]),
            for bi = 1:BlockNum
                disp(['Importing block: ',num2str(bi)]),
                % 读取并resize有效图片
                BlockFrameS = (bi-1)*BlockSz+1; % 当前block读取图片起始帧，从1开始编号
                BlockFrameE = min(bi*BlockSz,PeriodNumR(ci)); % 当前block读取图片结束帧，从1开始编号
                BlockFrameNum = BlockFrameE-BlockFrameS+1; % 当前block总图片帧数
                IMseriesT = zeros([ceil([pixR,pixC]),BlockFrameNum,ChAnlsNum(ci)],'uint16');
                for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
                    disp(['Importing(laser): ',ChTag{ci,Chfi}]),
                    IdxCamera = find(CorRatio.Camera==CameraTag{ci}); % 相机编号，用于找对应的CorRatio.data
                    IdxLaser = find(CorRatio.Laser==ChTag{ci,Chfi}); % 相机编号，用于找对应的CorRatio.data
                    ImCorRatio = imresize(CorRatio.data(:,:,IdxCamera,IdxLaser),1/ImBinning,'Method','nearest');   % 对应相机和激光条件下，CorRatio.data
                    tic,
                    for Prdi = BlockFrameS:BlockFrameE
                        if mod(Prdi,1000) == 0
                            disp(['Importing: ',num2str(Prdi)]),
                        end
                        %                 FileName = targetsName{(Prdi-1)*ChNum+ChAnalyseIdx(Chfi),ci};
                        FileName = targetsName{(PeriodNumRange(Prdi)-1)*ChAnlsNum(ci)+Chfi,ci};
                        % 读取原始图片
                        ImgRaw = double(imread(FileName)-RefDark(:,:,ci)); % 读取图片并减去相机offset
                        % 均匀性校正
                        ImgCor = ImgRaw.*ImCorRatio; % uniform/flatted correction
                        %                 figure,
                        %                 PixMin = min(ImgRaw,[],'all');
                        %                 PixMax = max(ImgRaw,[],'all');
                        %                 subplot(1,3,1);imshow(ImgRaw,[PixMin PixMax]),title('原始图像'); %显示原始图像
                        %                 subplot(1,3,2);imshow(ImCorRatio,[min(ImCorRatio(:)) max(ImCorRatio(:))]),title('校正系数'); %显示超像素分割图像
                        %                 subplot(1,3,3);imshow(ImgCor,[PixMin PixMax]),title('校正图像'); %显示超像素分割图像
                        %                 sgtitle(strrep(FileName,'_','\_'));
                        IMseriesT(:,:,Prdi-BlockFrameS+1,Chfi) = ImgCor;
                    end
                    toc,
                    % 平移校正
                    if MotionCorFlag
                        tic,
                        IMseriesTs = squeeze(IMseriesT(:,:,:,Chfi));
                        disp('Starting motion correction...');
                        if size(shiftsRcd,1) == PeriodNumR(ci) % 若已有完整shiftsRcd数据，则直接使用
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
            
            
            for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
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
            [optimizer, metric] = imregconfig('multimodal');%参数modality指定fixedIM image, movingIM image之间的关系，有两种选择‘monomodal’, 'multimodal'两种，
            optimizer.InitialRadius = optimizer.InitialRadius/3.5;%改变优化器的步长已达到对更加精细的变换
            % optimizer.MaximumIteCorns = 300;% 改变最大迭代次数
        end
        
        % 基于参考平均图片计算配准的转换矩阵
        if RegisterFlag
            fixedIM = IMseriesAve(:,:,fixCh,1);
            RfixedIM = imref2d(size(fixedIM));%imregtform把变化矩阵输出；然后用imref2d限制变换后的图像与参考图像有相同的坐标分布
            figure,imshow(fixedIM,[]);
            title(['fixedIM from Camera ',CameraTag{fixCh}]);
            for ci = 1:CameraNum % 同一个session根据每个相机记录的数据，CameraNum为相机数目
                if ci~=fixCh
                    movingIM = IMseriesAve(:,:,ci,1);
                    figure,imshow(movingIM,[]);
                    title(['movingIM of Camera ',CameraTag{ci}]);
                    figure(100),imshowpair(movingIM*3,fixedIM,'falsecolor'),title('Before registration');
                    tformSimilarity = imregtform(movingIM,fixedIM,'similarity',optimizer,metric);%用similarity的变换方式做初始配准，还可以用rigid，transform的方式
                    movingIMReg = imwarp(movingIM,tformSimilarity,'OutputView',RfixedIM);%imwarp函数执行几何变换，依据则是tformSimilarity的变换矩阵。
                    figure(101),imshowpair(movingIMReg,fixedIM,'falsecolor'),title('After registration');
                    %   saveas(gcf,[savefolder1,'\RegistCorn of ',imname(1:12),'.jpg']);
                    RegisterMat{ci,1} = tformSimilarity;
                    RefDarkRegister(:,:,ci) = imwarp(RefDark(:,:,ci),tformSimilarity,'OutputView',RfixedIM); %imwarp函数执行几何变换，依据则是tformSimilarity的变换矩阵。
                    
                    tic
                    IMseries(:,:,:,:,ci) = imwarp(IMseries(:,:,:,:,ci),tformSimilarity,'OutputView',RfixedIM); %imwarp函数执行几何变换，依据则是tformSimilarity的变换矩阵。
                    toc,
                    for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
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
        
        %% 保存IMseries
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
                WholeCortexMask(pixCi,pixRi) = nhoods(1,1); % 当Mask图片中单一像素点的值不同于八邻域像素点任一值时，用邻域像素点填充该点
            end
        end
    end
    figure,imshow(WholeCortexMask,[]);
    title('Mask imresized and filled');
    FileName = [MaskName(1:end-4),'_reScale.tif'];% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
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
            mask = mask.*WholeCortexMask;  % 非background ROI需要剔除血管
            mask0 = mask0.*WholeCortexMask0;  % 非background ROI需要剔除血管
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
    %% basal图片的帧数
    switch BasalFlagOfSession % 分析整个session所用的baseline图片生成方式：1-指定帧；2-根据刺激mark确定；3.整个session的每个像素点的X%低的数值；
        case 1 % 1-指定帧
            IdxBaseStart = IdxBaseList(fi,1); % 整个session的baseline开始帧数
            IdxBaseEnd = IdxBaseList(fi,2); % 整个session的baseline结束帧数
            if ~isnan(IdxBaseStart) & ~isnan(IdxBaseEnd) % 如果baseline未指定，则使用之前的IMseriesBasalIM
                IdxStart = max([IdxBaseStart,PeriodNumRange(1)])-PeriodNumRange(1)+1;
                IdxEnd = min([IdxBaseEnd,PeriodNumRange(end)])-PeriodNumRange(1)+1;
                if IdxBaseEnd==0
                    IdxEnd = ImSz3;
                end
            else
                IdxStart = IdxStart_bkp{fi-1}; % 使用其它session的baseline图片
                IdxEnd = IdxEnd_bkp{fi-1}; % 使用其它session的baseline图片
            end
            
        case 2  % 2-根据刺激mark确定，第一个刺激前的时间定为整个session的baseline结束时刻
            Baseline0Start = 0; % 整个session的baseline开始时刻
            Baseline0End = Stim.times(1,1); % 第一个刺激前的时间定为整个session的baseline结束时刻
            IdxStart = ceil(find(Imaging.times>=Baseline0Start,1)/ChNum/2);
            IdxEnd = floor(find(Imaging.times<=Baseline0End,1,'last')/ChNum/2);
            
        case 3  % 3-整个session的每个像素点的X%低的数值
            IdxBaseStart = IdxBaseList(fi,1); % 整个session的baseline开始帧数
            IdxBaseEnd = IdxBaseList(fi,2); % 整个session的baseline结束帧数
            IdxStart = max([IdxBaseStart,PeriodNumRange(1)])-PeriodNumRange(1)+1;
            IdxEnd = min([IdxBaseEnd,PeriodNumRange(end)])-PeriodNumRange(1)+1;
            
        case 4 % 4-根据刺激mark确定，session及所有trial的baseline并集
            IdxStart = cat(1,IdxBaseList(fi,1),IdxStartAll);
            IdxEnd = cat(1,IdxStartAll(1)-IdxCorLen,IdxStimAll-IdxCorLen);
            
        case 5 % 5-根据其它方式标注的baseline帧数确定，如睡眠的REM阶段,待改
            % Analysis by AccuSleep (1 = REM sleep, 2 = wakefulness, 3 = NREM sleep, 4 = undefined)
            StimRcdPath1 = fullfile(StimRcdPath,'SleepState');
            FileName = fullfile(StimRcdPath1,'SleepStateLabel.mat');
            load(FileName);
            
            %     labelsLen = size(labels,1);
            %     labelsTimeM = size(labels,1)*EpochTime;
            
            ImStamp.Flag = zeros(PeriodNumR(ci),1)*nan;
            ImStamp.Note = ImStamp.Flag;
            StimModeCount = zeros(max(labels),1);
            for Prdi = 1:PeriodNumR(ci) % 对每帧图像进行处理
                try
                    StimMode = labels(ceil(ImagingTimeM(Prdi)/EpochTime)); % 若最后一帧图像无对应时刻的睡眠状态划分，则跳过报错，使用前一个状态
                end
                ImStamp.Flag(Prdi,1) = StimMode;
                StimModeCount(StimMode,1) = StimModeCount(StimMode,1)+1;
                ImStamp.Note(Prdi,1) = StimModeCount(StimMode,1); % 该类型事件序号，即属于该事件/刺激的第几个trial
                ImStamp.Note(Prdi,2) = Prdi; % 总trial数编号
            end
            IdxStart = NaN;
            IdxEnd = NaN;
            FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_ImStamp.mat']);% 保存线性fit参数k，斜率，参数b，截距
            save(FileName,'ImStamp');
    end
    
    if BasalFlagOfSession == 5
        IdxBaseFrame = find(ImStamp.Flag==1); % 对睡眠觉醒数据，以REM对应时刻为baseline
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
        %% 减背景
        BackgroundRcd = zeros(CameraNum,ChAnlsNum(1),PeriodNumR(1),'uint16');
        ChNo = 0;
        figure,title('BackgroundRcd');hold on;
        for ci = 1:CameraNum % 同一个session根据每个相机记录的数据分别处理，CameraNum为相机数目
            for Chfi = 1:ChAnlsNum(ci) %  同一相机不同激发光
                IMseriesT = IMseries(:,:,:,Chfi,ci);
                IMseriesT = reshape(IMseriesT,[],PeriodNumR(ci));
                %     BackgroundRcd(ci,Chfi,:) = uint16(min(IMseriesT(mask_id{1},:),[],1)*BgToAutoF(ci,Chfi));
                BackgroundRcd(ci,Chfi,:) = uint16(mean(IMseriesT(mask_id{1},:),1));
                IMseriesT = IMseries(:,:,:,Chfi,ci)-BackgroundRcd(ci,Chfi,:);
                IMseries(:,:,:,Chfi,ci) = IMseriesT.*uint16(WholeCortexMask);  % 血管部分的像素点值设为0
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
        %% 基线校正（用于校正photobleach等基线的持续下降）
        if BaselineCorFlag
            % 调试发现每个像素点逐一fit,预计用时0.898313*38183/3600=9.5279h
            x = 1:ImSz3;
            BaseRange = 1:round(length(x)*0.05); % 使用前5%时间段的fit出的baseline作为baseline的预期真实值
            fdir = dir([MfileDir '\Funs\createFits.m']);
            if ~isempty(fdir)
                disp('createFits.m file exists');
            else
                %         yVarNameList = cell(CameraNum,ChAnlsNum(1));
                if BaselineCorFlag == 1
                    for ci = 1:CameraNum
                        for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
                            IMseriesT = IMseries(:,:,:,Chfi,ci);
                            IMseriesT = reshape(IMseriesT,[],size(IMseriesT,3));
                            y = mean(IMseriesT(mask_id{2},:),1);
                            yVarName = ['y',CameraTag{ci},ChTag{ci,Chfi}];
                            %                 yVarNameList{ci,Chfi} = yVarName;
                            eval([yVarName,'=y;']);
                            eval(['cftool(x,',yVarName,');']);
                        end
                    end
                    % 手动保存Matlab自动生成的createFits.m文件至mfile\Funs
                    helpdlg('Please generate createFits.m!');
                    pause;
                    uiwait(msgbox({'Please use sftool generate a function and saveas "createFits"';'Make sure ';'Finished?'},'Warning','modal'));
                end
            end
            
            % 对每个像素点进行baseline校正
            %           IMseries0 = IMseries;
            tic,
            disp('Baseline correction starting...');
            BlockParNum = 10; % 并行计算组数，分组以减小每个parfor占用的运行内存
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
                        %                          for PixjIdx = 1:length(PixjList)
                        %                                                 Pixi = 121;
                        %                                                 Pixj = 210;
                        %                         PixjIdx = 18;
                        Pixj = PixjList(PixjIdx);
                        if WholeCortexMask(Pixi,Pixj)
                            %                             tic,
                            disp(['Pix row:',num2str(Pixi),', Pix column:',num2str(Pixj)]);
                            y = zeros(4,PeriodNumR(1));
                            for ci = 1:CameraNum
                                for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
                                    Idx = (ci-1)*CameraNum+Chfi;
                                    y(Idx,:) = double(squeeze(IMseriesT(PixjIdx,:,Chfi,ci)));
                                end
                            end
                            if BaselineCorFlag == 1
                                [fitresult, gof] = createFits(x,  y(1,:),  y(2,:),  y(3,:), y(4,:));
                            end
                            yCorPix = zeros(PeriodNumR(1),ChAnlsNum(ci),CameraNum);
                            for ci = 1:CameraNum
                                for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
                                    Idx = (ci-1)*CameraNum+Chfi;
                                    %                               disp(Idx);
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
                                    %                               yCor = uint16(y(Idx,:)-yFit+mean(yFit));
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
                % RAM占用过高时清理并行池
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
            
            % 保存IMseries
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
    % 暂时生成basal图片用于后续计算
    IMseriesBasalIM = zeros(ImSz1,ImSz2,ChAnlsNum(1),CameraNum,'uint16');
    LowThd = 0.1; % 筛选每个像素点是否用于计算baseline的阈值
    for ci = 1:CameraNum
        for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
            [IMseriesBasalIM(:,:,Chfi,ci),IMseriesBasalIM_Thd(:,:,Chfi,ci)] = F_BasalIM_Ver1(IMseries(:,:,:,Chfi,ci),BasalFlagOfSession,LowThd,IdxBaseFrame);
            %         figure,imshow(IMseriesBasalIM(:,:,Chfi,ci),[]);
            % title(['Basal image: ',num2str(ci),', ',num2str(Chfi)]);
        end
    end
    %% 串光校正
    if SpectraUnmixFlag & CameraNum>1
        tic,
        disp('Spectra Unmixing starting...');
        %         % 用于生成二元一次方程组的解析式
        %                 syms realG realR Crg Cgr measuredG measuredR
        %                 eq1 = realG+realR*Crg-measuredG;
        %                 eq2 = realR+realG*Cgr-measuredR;
        %                 [realG,realR] = solve(eq1,eq2);
        %                 eval(realG)
        for Chfi = 1:ChAnlsNum(ci) %  同一相机不同激发光
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
    
    
    % 生成basal图片和AveSubBG图片
    IMseriesBasalIM = zeros(ImSz1,ImSz2,ChAnlsNum(1),CameraNum,'uint16');
    for ci = 1:CameraNum
        for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
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
    for ci = 1:CameraNum % 同一个session根据每个相机记录的数据分别处理，CameraNum为相机数目
        %% Artifacts correction
        tic,
        switch CorrectMethodFlag
            case 1
                % GB or RB Cor images
                disp('Artifacts correction by ratio image...');
                Chfi = 1; % 成像通道
                ChNormIdx = 2; % 405校正通道
                %             IMseriesMCor(:,:,:,Chfi,:) = double(IMseriesM(:,:,:,Chfi,:))./double(IMseriesM(:,:,:,ChNormIdx,:));
                IMseriesCor = single(IMseries(:,:,:,Chfi,ci))./single(IMseries(:,:,:,ChNormIdx,ci));
            case 2
                SaveFolder3 = [SaveFolder1,'\Regression(405nm) evaluation'];
                mkdir(SaveFolder3);
                disp('Artifacts correction by fitted 405 channel...');
                
                %                 IdxCorLen = 2; % 为避免光激活时基线校正使用的数据受光激活激发光的影响，所选数据点尽量前移，前移IdxCorLen个数据点
                %                 IdxStart = cat(1,IdxBaseList(fi,1),IdxStartAll);
                %                 IdxEnd = cat(1,IdxStartAll(1)-IdxCorLen,IdxStimAll-IdxCorLen);
                %                 IdxStart = IdxStart_bkp{fi};
                %                 IdxEnd = IdxEnd_bkp{fi};
                if BasalFlagOfSession == 5
                    IdxFit = find(ImStamp.Flag==2); % 对睡眠觉醒数据，以Wake对应时刻为fit的数据点
                else
                    IdxFit = IdxBaseFrame; % 基于baseline图片进行校正
                    %                     IdxFit0 = IdxBaseFrame; % 基于baseline图片进行校正
                    %                     IdxFit0 = []; % 原始记录数据中的图片编号
                    %                     for Idxi = 1:length(IdxStart)
                    %                         %                         IdxFit0 = cat(2,IdxFit0,IdxStart(Idxi):IdxEnd(Idxi));
                    %                         IdxFit0 = cat(2,IdxFit0,IdxStart(Idxi):IdxEnd(Idxi));
                    %                     end
                end
                %                 IdxFit = intersect(IdxFit0,PeriodNumRange)-PeriodNumRange(1)+1; % 在读取进来的IMseries中的图片编号
                
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
                    %                     img = squeeze(IMseriesT(:,:,IdxFit(i),1));
                    %                     IMseriesT(:,:,i,1) = imgaussfilt(squeeze(IMseriesT(:,:,IdxFit(i),1)),FilterSigma);
                    %                     img = squeeze(IMseriesT(:,:,IdxFit(i),2));
                    %                     IMseriesT(:,:,i,2) = imgaussfilt(squeeze(IMseriesT(:,:,IdxFit(i),2)),FilterSigma);
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
                            if mod(Pixi,50)==0 & mod(Pixj,50)==0 % 抽查像素点fit效果
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
                
                Chfi = 1; % 成像通道
                ChNormIdx =2; % 405校正通道
                IMseriesCor = single(IMseries(:,:,:,Chfi,ci))-thetaIMk.*single(IMseries(:,:,:,ChNormIdx,ci))-thetaIMb;
                IMseriesChNormBasalIM = IMseriesBasalIM(:,:,ChNormIdx,ci);
                IMseriesChNormBasalIMfit = thetaIMk.*single(IMseriesChNormBasalIM)+thetaIMb;
                figure,
                subplot(1,2,1),imshow(IMseriesChNormBasalIM,[]);title('Ex 405 nm');
                subplot(1,2,2),imshow(IMseriesChNormBasalIMfit,[]);title('Ex 405 nm fitted');
                FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},'_Ex405&fit images','.fig']);% 线性fit参数b，截距
                saveas(gcf,FileName);
                IMseriesCor = IMseriesCor+IMseriesChNormBasalIMfit; % 为避免F0过小，加上fit后的405数值
                clear IMseriesT;
        end
        toc,
        %% Analysis of response across the whole session
        if WholeSessionRspFlag
            disp('Analyzing response of the whole session...'),
            disp(['Generate basal images by method ',num2str(BasalFlagOfSession)]);
            % 输出Correction后的basal图片
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
            
            % #########输出单通道的反应图片#########
            compression = 0;
            IMseriesRsp = zeros(ImSz1,ImSz2,ImSz3,ImSz4,'single');
            tic,
            for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
                IMseriesRsp(:,:,:,Chfi) = squeeze(single(IMseries(:,:,:,Chfi,ci)).*single(WholeCortexMask)./single(IMseriesBasalIM(:,:,Chfi,ci)))-1;
            end
            IdxBad = find(abs(IMseriesRsp)==Inf);
            if IdxBad
                IMseriesRsp(IdxBad) = 0; % 将无穷大的异常值置零，避免影响整体平均值
            end
            for Chfi = 1:ChAnlsNum(ci) %1:ChNum  %读取同一相机记录文件夹下的图片
                if SaveSingleChFlag(ci)
                    FileName = fullfile(SaveFolder1,[SessionChTitle{fi,ci},'_',ChTag{ci,Chfi},'_Rsp.tif']);
                    F_WriteBigTiff_Ver1(FileName,pixelsize,compression,IMseriesRsp(:,:,:,Chfi));
                    if BasalFlagOfSession == 5
                        % 输出不同睡眠状态的平均图片及mat文件
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
            % #########输出Cor的反应图片#########
            tic,
            IMseriesRspCor = IMseriesCor.*single(WholeCortexMask)./IMseriesCorBasalIM(:,:,ci)-1;
            IdxBad = find(abs(IMseriesRspCor)==Inf);
            if IdxBad
                IMseriesRsp(IMseriesRspCor) = 0; % 将无穷大的异常值置零，避免影响整体平均值
            end
            Chfi = 1;
            FileName = fullfile(SaveFolder1,[SessionTitle{fi,1},'_',CameraTag{ci},'_',ChTag{ci,Chfi},ChTag{ci,ChNormIdx},'cor_Rsp','.tif']); % num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
            F_WriteBigTiff_Ver1(FileName,pixelsize,compression,IMseriesRspCor);
            if BasalFlagOfSession == 5
                % 输出不同睡眠状态的平均图片及mat文件
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
            RspOfSession.ChS = cell(i,size(ROIs,2)-1); % 记录通道,cell(align种类数,ROI数目)
            RspOfSession.ChNorm = cell(i,size(ROIs,2)-1); % 405校正通道
            RspOfSession.Cor = cell(i,size(ROIs,2)-1); % 校正后通道
            RspOfSession.StartTimeInSec = ImagingTimeM(1);
            RspOfSession.TimeInSec = ImagingTimeM'-RspOfSession.StartTimeInSec;
            RspOfSession.TimeInMin = RspOfSession.TimeInSec/60;
            RspOfSession.TimeInHrs = RspOfSession.TimeInSec/3600;
            MFIOfSession.ChS = cell(i,size(ROIs,2)-1); % 记录通道的平均荧光亮度
            MFIOfSession.ChNorm = cell(i,size(ROIs,2)-1); % 校正通道平均荧光亮度
            MFIOfSession.Cor = cell(i,size(ROIs,2)-1); % 校正后的平均荧光亮度
            
            % ########记录通道########
            IMseriesRspChS = IMseriesRsp(:,:,:,1);
            IMseriesRspChS = reshape(IMseriesRspChS,[],PeriodNumR(ci));
            %             IMseriesRspChS = IMseriesRspChS';
            % ########405校正通道########
            IMseriesRspChNorm = IMseriesRsp(:,:,:,2);
            IMseriesRspChNorm = reshape(IMseriesRspChNorm,[],PeriodNumR(ci));
            %             IMseriesRspChNorm = IMseriesRspChNorm';
            % ########校正后通道########
            IMseriesRspCor = reshape(IMseriesRspCor,[],PeriodNumR(ci));
            %             IMseriesRspCor = IMseriesRspCor';
            if MFIrecordFlag  % 记录ROI随时间变化的平均荧光亮度
                % ########记录通道########
                IMseriesMFIChS = IMseries(:,:,:,1,ci);
                IMseriesMFIChS = reshape(IMseriesMFIChS,[],PeriodNumR(ci));
                %                 IMseriesMFIChS = IMseriesMFIChS';
                % ########405校正通道########
                IMseriesMFIChNorm = IMseries(:,:,:,2,ci);
                IMseriesMFIChNorm = reshape(IMseriesMFIChNorm,[],PeriodNumR(ci));
                %                 IMseriesMFIChNorm = IMseriesMFIChNorm';
                % ########校正后通道########
                IMseriesCorMFI = reshape(IMseriesCor,[],PeriodNumR(ci));
                %                 IMseriesCorMFI = IMseriesCorMFI';
            end
            
            for ri = 2:size(ROIs,2)
                TempValue = mean(IMseriesRspChS(mask_id{ri},:),1,'omitnan');
                RspOfSession.ChS{i,ri-1} = cat(2,RspOfSession.ChS{i,ri-1},TempValue'); % 记录通道的反应
                TempValue = mean(IMseriesRspChNorm(mask_id{ri},:),1,'omitnan');
                RspOfSession.ChNorm{i,ri-1} = cat(2,RspOfSession.ChNorm{i,ri-1},TempValue'); % 校正通道反应
                TempValue = mean(IMseriesRspCor(mask_id{ri},:),'omitnan');
                RspOfSession.Cor{i,ri-1} = cat(2,RspOfSession.Cor{i,ri-1},TempValue'); % 比值的反应
                if MFIrecordFlag
                    TempValue = mean(IMseriesMFIChS(mask_id{ri},:),1,'omitnan');
                    MFIOfSession.ChS{i,ri-1} = cat(2,MFIOfSession.ChS{i,ri-1},TempValue'); % 记录通道的平均荧光亮度
                    TempValue = mean(IMseriesMFIChNorm(mask_id{ri},:),1,'omitnan');
                    MFIOfSession.ChNorm{i,ri-1} = cat(2,MFIOfSession.ChNorm{i,ri-1},TempValue'); % 校正通道平均荧光亮度
                    TempValue = mean(IMseriesCorMFI(mask_id{ri},:),1,'omitnan');
                    MFIOfSession.Cor{i,ri-1} = cat(2,MFIOfSession.Cor{i,ri-1},TempValue'); % 比值的平均荧光亮度
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
        if StimRspFlag % 根据刺激信息分析整理图片
            tic,
            StimTag = StimTagList{fi};   % 刺激名称
            RspOfSingleTrial.ChS = cell(StimModeNum,size(ROIs,2)-1); % 记录通道
            RspOfSingleTrial.ChNorm = cell(StimModeNum,size(ROIs,2)-1); % 405校正通道
            RspOfSingleTrial.Cor = cell(StimModeNum,size(ROIs,2)-1); % 比值通道
            Chfi = 1; % 成像通道
            ChNormIdx = 2; % 405校正通道
            for i = 1:StimModeNum
                if StimModeCount(i)>0
                    IMseriesTrialRspChS{i} = zeros([ImSz1,ImSz2,TrialFrames,StimModeCount(i)],'single');
                    IMseriesTrialRspChNorm(i) = IMseriesTrialRspChS(i);
                    IMseriesTrialRspCor(i) = IMseriesTrialRspChS(i);
                    IdxStim = find(ImStamp.Flag(:,1) == i); % 距离刺激编号为i的最近一帧图片的编号
                    IdxTrial = ImStamp.Note(IdxStim,2); % 每一个trial对应总trial数的第IdxTrial个trial
                    IdxStim = IdxStim-PeriodNumRange(1)+1;
                    IdxStart = IdxStartAll(IdxTrial);
                    IdxEnd = IdxEndAll(IdxTrial);
                    for ti = 1:StimModeCount(i) % 该类型trial总数
                        IMseriesTrialChS = single(IMseries(:,:,IdxStart(ti):IdxEnd(ti),Chfi,ci)); % 感兴趣通道图片提取
                        %                         IMbasal = mean(IMseriesTrialChS(:,:,1:IdxStim(ti)-IdxStart(ti),Chfi),3);
                        IMbasal = mean(IMseriesTrialChS(:,:,1:IdxStim(ti)-IdxStart(ti)-IdxCorLen,Chfi),3);
                        RspChS  = IMseriesTrialChS./IMbasal-1;
                        IMseriesTrialRspChS{i}(:,:,:,ti) = RspChS;
                        IMseriesTrialChNorm = single(IMseries(:,:,IdxStart(ti):IdxEnd(ti),ChNormIdx,ci)); % 校正通道图片提取
                        IMbasal = mean(IMseriesTrialChNorm(:,:,1:IdxStim(ti)-IdxStart(ti)-IdxCorLen,Chfi),3);
                        RspChNorm = IMseriesTrialChNorm./IMbasal-1;
                        IMseriesTrialRspChNorm{i}(:,:,:,ti) = RspChNorm;
                        % 在未计算整个session的反应的时候，为减少内存占用而使用以下代码
                        switch CorrectMethodFlag
                            case 1
                                IMseriesTrialCor = IMseriesTrialChS./IMseriesTrialChNorm;
                            case 2
                                IMbasalfit = thetaIMk.*IMbasal+thetaIMb;
                                IMseriesTrialCor = IMseriesTrialChS-thetaIMk.*IMseriesTrialChNorm-thetaIMb+IMbasalfit;
                        end
                        %                         figure,imshow(IMbasal,[]);
                        %
                        %                                 RspCor = IMseriesTrialCor./IMbasal-1;figure,imshow(mean(IMseriesTrialCor,3),[]);
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
                        % 488nm or 561nm通道
                        Value = RspOfSingleTrial.ChS{i,ri-1}; % 反应值
                        ValueM = mean(Value,2); % 平均反应
                        ValueSEM = std(Value,0,2,'omitnan')/sqrt(size(Value,2)); % SEM
                        RspOfSingleTrial.ChS_M{i,ri-1} = cat(2,ValueM,ValueSEM);
                        % 405nm校正通道
                        Value = RspOfSingleTrial.ChNorm{i,ri-1}; % 反应值
                        ValueM = mean(Value,2); % 平均反应
                        ValueSEM = std(Value,0,2,'omitnan')/sqrt(size(Value,2)); % SEM
                        RspOfSingleTrial.ChNorm_M{i,ri-1} = cat(2,ValueM,ValueSEM);
                        % 488nm/405nm校正通道
                        Value = RspOfSingleTrial.Cor{i,ri-1}; % 反应值
                        ValueM = mean(Value,2); % 平均反应
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
            clear IMseriesTrialRspCor RspOfSingleTrial IMseriesTrialRspChNorm IMseriesTrialRspChS; % 释放内存
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
%% 符合特定名称的文件或者文件夹扫描函数scanDir
function TargetList = scanDir(CurrentPath,TargetLabel,TargetType,MatchPattern)
% 遍历CurrentPath下的文件夹及文件
% CurrentPath = mainPath;
files = dir(CurrentPath);
len = length(files);
TargetList = {};
index = 1;
for i = 1:len
    % 跳过.以及..文件夹
    if (strcmp(files(i).name, '.') == 1) ...
            || (strcmp(files(i).name, '..') == 1)
        continue;
    end
    % 递归调用函数，遍历当前目录下的文件夹(深度过深，可能会报错)
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
% 用于输出大tiff文件的函数，可以输出大于4GB的tiff文件
% Fei Deng,20210808
function F_WriteBigTiff_Ver1(FileName,pixelsize,compression,ImStack)
obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
[ImSz1 ImSz2 ImSz3] = size(ImStack);
for i = 1:ImSz3
    img = ImStack(:,:,i);
    img = img'; % 虽然现在的img和实际的方向一致，但是存储以后不一致，因而需要转置以后再保存
    WriteIMG(obj,img);
end
close(obj);
end
%% =======================================================================================
%% F_StimModeImport_Ver0,用于读取txt文本中的StimMode记录
%## Fei Deng,20210620,用于读取txt文本中的StimMode记录
function StimMode = F_StimModeImport_Ver0(fullFileName,SplitTag)
% SplitTag = 'StimMode: ';
fid = fopen(fullFileName,'rt'); % t是告诉fread是这里文本文件
% 统计txt行数
txtRow = 0;
while ~feof(fid)
    txtRow = txtRow+sum(fread(fid,10000,'*char')==char(10));
    % 一次性读取10000字符，计算其中的回车个数，其中10是回车的ASCII编码
    % '*char'表示每次读取一个字符，*表示输出也是字符
    % fread现在已经可以自动识别中文了，万一还是识别不了，请在fopen中指定文件编码格式，比如gbk
end
fclose(fid);
disp(['There are ',num2str(txtRow),' lines in txt file, importing data:']);
% 读取txt原始数据
recordV = cell(txtRow,2); % 预先分配存储空间，提高运算速度
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
    %  m = split(tline,'->'); % 分离每行的列
    m = split(tline,SplitTag); % 分离每行的列
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
%% ## Fei Deng,20220117,用于输入Dark images，生成offset图片
function RefDark = F_RefDarkImport_Ver0(TargetPathD,filetype,pixR,pixC,CameraNum,pixelSzBin1,ImBinning,compression)
disp('Looking for dark images in folder:');
disp(TargetPathD);
RefDark = uint16([]);
try
    DarkDir = fullfile(TargetPathD,'Dark');
    cd(DarkDir);
    DarkList = dir('*Dark*');
    if size(DarkList,1) == CameraNum
        for ci = 1:CameraNum  % 同一个session根据每个相机记录的数据分别处理，CameraNum为相机数目
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
%% ## Fei Deng,20220328,用于生成IMseries3DBasalIM图片
% Fei Deng,20220413,增加对输入数据类型的自动判断
%% =======================================================================================
function [BasalIM,BasalIM_Thd] = F_BasalIM_Ver1(IMseries3D,BasalFlagOfSession,LowThd,IdxBaseFrame)
className = class(IMseries3D); % 自动识别数据类型
[ImSz1,ImSz2,ImSz3] = size(IMseries3D);
BasalIM = zeros(ImSz1,ImSz2);
BasalIM_Thd = BasalIM;
% IMseries3DBasalIM = zeros(ImSz1,ImSz2),ImSz3,'uint16');
switch BasalFlagOfSession
    case 3
        tic,
        %             LowThd = 0.1; % 筛选每个像素点是否用于计算baseline的阈值
        for coli = 1:ImSz2
            %                                     coli = 150
            img = squeeze(IMseries3D(:,coli,:)); % 按列降维得到每一列像素点随时间变化的图像
            %                             figure,imshow(img,[]);
            imgSort = sort(img,2);
            %                             figure,imshow(imgSort,[]);
            BasalIM(:,coli) = mean(imgSort(:,1:round(LowThd*ImSz3)),2); % imgSort的每一行1:round(LowThd*ImSz3)小的数据点均值，每一列像素点的各自baseline值
            BasalIM_Thd(:,coli) = imgSort(:,round(LowThd*ImSz3)); % 每一列像素点的各自baseline阈值
        end
    otherwise
        BasalIM = squeeze(mean(IMseries3D(:,:,IdxBaseFrame),3));
end
BasalIM = eval([className,'(BasalIM)']); % 更改数据类型和输入数据类型一致
BasalIM_Thd = eval([className,'(BasalIM_Thd)']); % 更改数据类型和输入数据类型一致
% figure,imshow(BasalIM,[]);
% title('Basal image');
% close;
end

%   =======================================================================================
% 程序修改记录
% Fei Deng,20210517,for response analysis with stimuli
% Fei Deng,20210607,修复了WriteTiff输出tiff文件有时出错的问题，改用try语句反复调用函数分段写入
% Fei Deng,20210607,修正了IMseriesBasalM的图片帧数挑选
% Fei Deng,20210608,修正了原始数据路径TargetPath对应关系
% Fei Deng,20210620,适用于多种刺激模式
% Fei Deng,20210725,适用于单次刺激内部包含编码信息，根据trial时长进行判定
% Fei Deng,20210806,% 针对spike2记录的帧数比实际成像多1帧的debug，找出多出来的数据点并删除
% Fei Deng,20210808,% 改进后可以输出大于4GB的tiff文件
% Fei Deng,20210906,% 所有输出图片都添加scale => StimRspAnalysis_Ver5b_For1Ch
% Fei Deng,20210906,% 增加两个不同相机之间图片的配准 => StimRspAnalysis_Ver6_For2Ch
% Fei Deng,20210906,% 增加ROI绝对亮度MFI的记录
% Fei Deng,20211101,% 细节优化完善中,使用SleepWake数据调试
% Fei Deng,20211122,% 使用SleepWake数据完成调试
% Fei Deng,20211126,% 优化对Dark图片的读取，已有平均后的dark imagea便直接读取
% Fei Deng,20211126,% 增加对原始数据分析帧数范围输入
% Fei Deng,20211213,% 增加对每个session的baseline帧数范围输入
% Fei Deng,20211216,% 通过mfile文件名识别待处理原始数据文件夹名及处理结果文件夹名
% Fei Deng,20211227,% 增加对spike2记录时检测到的imaging通道异常高电平进行自动剔除
% Fei Deng,20220103,% 增加输出resize后的mask图片
%   =======================================================================================
