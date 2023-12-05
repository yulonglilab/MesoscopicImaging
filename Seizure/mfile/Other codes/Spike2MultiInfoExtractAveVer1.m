%% Spike2MultiInfoExtractAveVer1
% Running time:20221019
%   ======================================================================================
% Fei Deng,20221019,基于Spike2MultiInfoExtractVer6d修改，用于平均多种刺激信息，视频输出部分暂未调试
% =======================================================================================
%   数据说明：
%   时间单位s
% =======================================================================================
close all;
clearvars -except CropList,clc;
% clear,clc,
ROI_No = 1; % 分析的感兴趣的ROI
VideoFlag = 0; % 是否截取视频
LabelFlag = 1; % 是否添加文字标注
TimeFlag = 1; % 视频标注中时间显示模式：1- s, 2-MM:SS, 3-HH:MM:SS
FrameInterNum = 1; % 生成response video的取图片间隔，每FrameInterNum取一帧
FrameRate = 10; % 输出的response video的帧率，FrameRate帧/s
% VideoTimeList = [25;25;70]; % 手动填写输出视频时长（s）
FrameIntvlStep = 5; % 输出视频时间尺度上间隔多少帧取一张，用于减小输出文件尺寸
BrightIndex = [1.2,1]; % 视频图像亮度调整系数
TargetLabel = 'IMseriesOfSingleTrial.mat'; % 输入目标文件夹/文件特征字符串
%% imaging数据校正
% BadIdx = 1:52; % 需要修正的数据范围
% Imaging0 = Imaging;
% Imaging.times(BadIdx) = [];
% Imaging.level(BadIdx) = [];
% Imaging.length = length(Imaging.level);

%% 根据成像提取spike2记录信息
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
cd(ParentFolder);
% [file,path] = uigetfile('*RspOfSession.mat','请选择 RspOfSession.mat文件'); %文件路径
[file,path] = uigetfile('*IMseriesOfSingleTrial.mat','请选择 IMseriesOfSingleTrial.mat文件'); %文件路径
ChTagIdx = strfind(file,TargetLabel)-2;
SaveFolder = [path,'Synchronized data'];
savetime = datestr(now,'yyyymmddTHHMMSS');
SaveFolder = [SaveFolder,'_',savetime];
mkdir(SaveFolder);
TargetFileName = file;
TargetFileName(ChTagIdx) = '*';
cd(path);
targetPack = dir(TargetFileName); %待整理文件列表获取
ChNum = size(targetPack,1); % 待转换的帧数编号，自动计算
FileName = fullfile(path,file);
load(FileName);
ImagingTimeMs = RspOfSingleTrial.TimeInSec;
% load *Stim&Imaging Info.mat
file = dir('*Stim&Imaging Info.mat'); %文件路径
FileName = file.name;
% load(FileName);
load(FileName,'IdxStartAll','IdxEndAll','IdxStimAll','Image_f','ImagingTimeM','StimMode','StimModeNum','ImStamp');
BlockListS = IdxStartAll;
BlockListE = IdxEndAll;
BlockFrameNumList = cat(2,BlockListS,BlockListE); % 待截取block的帧数范围

VideoTimeList = (BlockFrameNumList(:,2)-BlockFrameNumList(:,1)+1)/FrameInterNum/FrameRate; % 自动计算输出视频时长（s）

CropList = [];
%% 提取每个block的数据
StimModeCountEach = zeros(StimModeNum,1);
for fi = 1:size(targetPack,1)
    FileName = targetPack(fi).name;
    ChTag = FileName(ChTagIdx);
    FilePathName = fullfile(path,FileName);
    load(FilePathName);
    for bi = 1:size(BlockFrameNumList,1)
        StimModeTemp = StimMode(bi);
        StimModeCountTemp = ImStamp.Note(IdxStimAll(bi));
        BlockFrameNumS = BlockFrameNumList(bi,1); % 待截取block的起始帧数
        BlockFrameNumE = BlockFrameNumList(bi,2); % 待截取block的终止帧数
        %         ImPeriodT = (ImagingTimeM(end)-ImagingTimeM(1))/(length(ImagingTimeM)-1);
        %         ImagingTimeMb = ImagingTimeM(BlockFrameNumS:BlockFrameNumE)-ImagingTimeM(BlockFrameNumS);
        ImagingTimeMb = ImagingTimeMs;
        %     TableSync = {'Time (s)','Time (min)','Time (h)'};
        %     ImagingTimeMb_hms = repmat(ImagingTimeMb,1,3)./[1,60,3600];
        %         TableSync(2:size(ImagingTimeMb,1)+1,:) = num2cell(ImagingTimeMb_hms);
        %     filename = fullfile(SaveFolder,['ImagingTime_im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'.csv']);
        %     writecell(TableSync,filename);
        RspOfSingleTrialBlock = cat(2,ImagingTimeMb,RspOfSingleTrial.Cor{StimModeTemp, ROI_No}(:,StimModeCountTemp),RspOfSingleTrial.ChS{StimModeTemp, ROI_No}(:,StimModeCountTemp),RspOfSingleTrial.ChNorm{StimModeTemp, ROI_No}(:,StimModeCountTemp));
        TableSync = {'Time (s)','Cor','ChS','ChNorm'};
        TableSync(2:size(ImagingTimeMb,1)+1,:) = num2cell(RspOfSingleTrialBlock);
        TrialTitleCh = FileName(1:strfind(FileName,'IMseriesOfSingleTrial.mat')-1);
        TrialTitle = TrialTitleCh(1:end-2);
        filename = fullfile(SaveFolder,[TrialTitleCh,'im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'.csv']);
        disp(['Saving file: ',filename]);
        writecell(TableSync,filename);
        disp('File saved!');
        if fi == 1
            %         RspOfSingleTrial.StartTimeInSec = 1.69; % !!!!!!!!!!!!!!!! 手动添加
%             BlockStartT = ImagingTimeM((BlockFrameNumS))+ImagingTimeM(1); % 所截取部分在整个session的起始时刻（s）
%             BlockEndT = ImagingTimeM((BlockFrameNumE))+ImagingTimeM(1); % 所截取部分在整个session的终止时刻（s）
            BlockStartT = ImagingTimeM((BlockFrameNumS)); % 所截取部分在整个session的起始时刻（s）,ImagingTimeM(1)即为spike2记录的真实时刻
            BlockEndT = ImagingTimeM((BlockFrameNumE)); % 所截取部分在整个session的终止时刻（s）,ImagingTimeM(1)即为spike2记录的真实时刻
            EEGb = EEG;
            EEGb.values = EEG.values(max([ceil(BlockStartT/EEG.interval),1]):min([ceil(BlockEndT/EEG.interval),EEG.length]));
            EEGb.length = length(EEGb.values);
            AllEEGraw{StimModeTemp,StimModeCountTemp} = EEGb; % 存储每个trial的数据至cell类型中，行对应刺激模式编号，列对应重复数编号
            AllData.EEG{StimModeTemp,1}(1:EEGb.length,StimModeCountTemp) = EEGb.values; % 存储每个trial的数据至cell类型中，行对应刺激模式编号，cell内部矩阵行对应时间维度的数据点，列对应重复数编号
            
            EMGb = EMG;
            EMGb.values = EMG.values(max([ceil(BlockStartT/EMG.interval),1]):min([ceil(BlockEndT/EMG.interval),EMG.length]));
            EMGb.length = length(EMGb.values);
            AllEMGraw{StimModeTemp,StimModeCountTemp} = EMGb; % 存储每个trial的数据至cell类型中，行对应刺激模式编号，列对应重复数编号
            AllData.EMG{StimModeTemp,1}(1:EMGb.length,StimModeCountTemp) = EMGb.values; % 存储每个trial的数据至cell类型中，行对应刺激模式编号，cell内部矩阵行对应时间维度的数据点，列对应重复数编号
            
            Runb = Run;
            Runb.values = Run.values(max([ceil(BlockStartT/Run.interval),1]):min([ceil(BlockEndT/Run.interval),Run.length]));
            Runb.length = length(Runb.values);
            % run speed
            MAXDACCNTS = 4095; % maximum dac value
            ZeroDACVOLTS = 1.5; % 速度为0对应的输出电压
            MAXSPEED = 1000; % maximum speed for dac out (mm/sec),此处单位为mm
            maxDACval = 3; % DAC ouput voltage at maximum speed
            ZeroDACval = ZeroDACVOLTS * MAXDACCNTS / 3.3;
            Velocity = Runb.values;
            dacval = Velocity/3.3*MAXDACCNTS;
            % Velocity = (dacval - ZeroDACval)/ZeroDACval*MAXSPEED; % 计算绝对速度，单位cm/s
            Velocity = (dacval - ZeroDACval)/ZeroDACval*MAXSPEED/10; % 计算绝对速度，单位cm/s，即速度=(实际电压-速度零对应电压)/速度零对应电压*100 cm/s = (实际电压/1.5-1)*100cm/s;
            Runb.valuesInCM = Velocity;
            AllRunraw{StimModeTemp,StimModeCountTemp} = Runb; % 存储每个trial的数据至cell类型中，行对应刺激模式编号，列对应重复数编号
            AllData.Run{StimModeTemp,1}(1:Runb.length,StimModeCountTemp) = Runb.valuesInCM; % 存储每个trial的数据至cell类型中，行对应刺激模式编号，cell内部矩阵行对应时间维度的数据点，列对应重复数编号
            
            Timeb = Run.interval*(1:Runb.length)';
            TableSync = {'Time (s)','Velocity (cm/s)'};
            % Timeb = Timeb/60;
            % TableSync = {'Time (min)','Velocity (cm/s)'};
            Rundata = cat(2,Timeb,Velocity);
            TableSync(2:Runb.length+1,:) = num2cell(Rundata);
            filename = fullfile(SaveFolder,[TrialTitle,'im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'_Velocity.csv']);
            disp(['Saving file: ',filename]);
            writecell(TableSync,filename);
            disp('File saved!');
            
            Timeb = EEGb.interval*(1:EEGb.length)';
            TableSync = {'Time (s)','EEG (uV)','EMG (uV)'};
            % Timeb = Timeb/60;
            % TableSync = {'Time (min)','EEG (μV)','EMG (μV)'};
            EEGEMGdata = cat(2,Timeb,EEGb.values,EMGb.values);
            TableSync(2:EEGb.length+1,:) = num2cell(EEGEMGdata);
            filename = fullfile(SaveFolder,[TrialTitle,'im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'_EEGEMG.csv']);
            writecell(TableSync,filename);
            BinStep = 10;
            BinLen = floor(size(EEGEMGdata,1)/BinStep);
            EEGEMGdataBin = zeros(BinLen,3);
            for i = 1:BinLen
                IdxStart = (i-1)*BinStep+1;
                IdxEnd = i*BinStep;
                EEGEMGdataBin(i,:) = mean(EEGEMGdata(IdxStart:IdxEnd,:),1);
            end
            TableSync = {'Time (s)','EEG (uV)','EMG (uV)'};
            TableSync(2:size(EEGEMGdataBin,1)+1,:) = num2cell(EEGEMGdataBin);
            filename = fullfile(SaveFolder,[TrialTitle,'im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'_EEGEMG_Bin Step=',num2str(BinStep),'.csv']);
            disp(['Saving file: ',filename]);
            writecell(TableSync,filename);
            disp('File saved!');
            %% 截取行为视频
            if VideoFlag == 1  % 是否截取视频
                cd(StimRcdFolder);
                CropFlag = 1;
                filetype = 'mp4';
                Videofile = dir(['*.',filetype]); %文件路径
                frameStartT = BlockStartT;
                frameEndT = BlockEndT;
                for i = 1:size(Videofile,1)
                    VideoT = VideoTimeList(bi);
                    VideoName = Videofile(i).name;
                    disp(['Editing video: ',VideoName]);
                    VideoR = VideoReader(VideoName);
                    % while hasFrame(v)
                    %     frame = readFrame(v);
                    % end
                    
                    if frameStartT == 0
                        frameStart = 1;
                    else
                        frameStart = round(frameStartT*VideoR.FrameRate);
                    end
                    if frameEndT == 0
                        frameEnd = VideoR.NumFrames;
                    else
                        frameEnd = round(frameEndT*VideoR.FrameRate);
                    end
                    
                    % 读取视频中待展示图片
                    VideoBlock1st = read(VideoR,frameStart);
                    [Sz1 Sz2 Sz3] = size(VideoBlock1st);
                    FrameRange = frameStart:FrameIntvlStep:frameEnd;
                    Sz4 = length(FrameRange);
                    VideoBlock = repmat(VideoBlock1st,[1 1 1 Sz4]);
                    disp('Importing images from video...');
                    tic,
                    for j = 1:Sz4
                        if mod(j,100) == 0
                            disp([num2str(j),'/',num2str(Sz4)]);
                        end
                        VideoBlock(:,:,:,j) = read(VideoR,FrameRange(j));
                    end
                    toc
                    
                    figure(101),imshow(VideoBlock(:,:,:,1)*BrightIndex(i),[]);
                    %                     figure,imshow(VideoBlock(:,:,:,1)*1.2,[]);
                    title('Raw');
                    %                figure(102),imshow(VideoBlock(:,:,:,1)*1.3,[]);
                    %             title('Brightness adjust');
                    if CropFlag
                        if size(CropList,1)<i
                            [CropX,CropY] = ginput(2); % 在图中标记矩形截取区域的左上角和右下角顶点
                            CropX = [max(uint16(CropX(1)),1),min(uint16(CropX(2)),Sz2)]; % 避免选到边缘报错
                            CropY = [max(uint16(CropY(1)),1),min(uint16(CropY(2)),Sz1)];
                            
                            CropList(i,1:4) = [CropX,CropY]; % 暂存每个行为学相机的ROI
                        else
                            CropX = CropList(i,1:2); % 读取已经暂存的每个行为学相机的ROI
                            CropY = CropList(i,3:4);
                        end
                        figure(102),imshow(VideoBlock(CropY(1):CropY(2),CropX(1):CropX(2),:,1),[]);
                        title('Cropped');
                        framesCrop = VideoBlock(CropY(1):CropY(2),CropX(1):CropX(2),:,:);
                        [Sz1 Sz2 Sz3 Sz4] = size(framesCrop);
                    else
                        framesCrop = VideoBlock;
                    end
                    %             framesCrop = framesCrop*BCindex(i);
                    %             figure(103),
                    %             for j = 1:Sz4
                    %                 disp(j);
                    %                 imshow(framesCrop(:,:,1,j),[]);
                    %             end
                    %% 添加标注
                    framesCropNoted = framesCrop;
                    if LabelFlag
                        PeriodOfVideoW = FrameIntvlStep/VideoR.FrameRate;
                        posX = round(0.025*Sz2);
                        posY = round(0.015*Sz1);
                        position = [posX,posY;posX,posY;posX,posY];
                        box_color = {'red','green','yellow'};
                        for j = 1:Sz4
                            FrameCurTime = PeriodOfVideoW*j;
                            switch TimeFlag
                                case 1
                                    text_str = [num2str(round(FrameCurTime)),' s'];
                                case 2
                                    text_str = datestr(FrameCurTime/3600/24,'MM:SS');
                                case 3
                                    text_str = datestr(FrameCurTime/3600/24,'HH:MM:SS');
                            end
                            %     box_color = {'red','green','yellow'};
                            framesCropNoted(:,:,:,j) = insertText(framesCrop(:,:,:,j)*BrightIndex(i),position,text_str,...
                                'FontSize',30,'BoxColor',box_color,'BoxOpacity',0,'TextColor','white');
                            %     figure,imshow(framesCropNoted(:,:,:,j));
                        end
                    end
                    %% 输出视频
                    % VideoNameW = [VideoName(1:strfind(VideoName,filetype)-2),'_frame',num2str(frameStart),'-',num2str(frameEnd),'-',num2str(FrameIntvlStep),'_',num2str(VideoT),'s.avi'];
                    %             VideoNameW = fullfile(SaveFolder,[VideoName(1:strfind(VideoName,'mp4')-2),'_T',num2str(frameStartT),'-',num2str(frameEndT),'s-',num2str(FrameIntvlStep),'_',num2str(VideoT),'s.avi']);
                    VideoNameW = fullfile(SaveFolder,[VideoName(1:strfind(VideoName,'mp4')-2),'_im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'_T',num2str(frameStartT),'-',num2str(frameEndT),'s-',num2str(FrameIntvlStep),'_',num2str(VideoT),'s.avi']);
                    VideoW = VideoWriter(VideoNameW,'Motion JPEG AVI');
                    VideoW = VideoWriter(VideoNameW);
                    % VideoW = VideoWriter(VideoNameW,'Uncompressed AVI');
                    % VideoW.Quality  = 100;
                    VideoW.FrameRate = round(Sz4/VideoT);
                    open(VideoW);
                    disp(['Saving file: ',filename]);
                    writeVideo(VideoW,framesCropNoted);
                    disp('File saved!');
                    % writeVideo(VideoW,frames);
                    close(VideoW);
                    close(101);
                    close(102);
                    clear framesCrop VideoBlock frameStart frameEnd;
                    disp('Finished!');
                end
            end
        end
    end
end
%% 整理汇总后的数据
for i = 1:StimModeNum
    TempData = AllData.EEG{i,ROI_No};
    TrialNum = size(TempData,2);
    % 将末尾的连续0置为NaN
    for j = 1:TrialNum
        IdsGood(j) = find(TempData(:,j),1,'last');
        TempData(IdsGood(j)+1:end,j) = NaN;
    end
    AllData.EEG_M{i,1}(:,1) = mean(TempData,2,'omitnan'); % 平均反应
    AllData.EEG_M{i,1}(:,2) = std(TempData,0,2,'omitnan')/sqrt(size(TempData,2)); % SEM
    %     AllData.EEG_SEM{i,1} = std(TempData,0,2,'omitnan')/sqrt(size(TempData,2)); % SEM
    AllData.EEG_Time{i,1} = AllEEGraw{1, 1}.interval*(1:size(TempData,1))'; % 时间（s）
    AllData.EEG{i,ROI_No} = TempData;
    
    TempData = AllData.EMG{i,ROI_No};
    TrialNum = size(TempData,2);
    % 将末尾的连续0置为NaN
    for j = 1:TrialNum
        IdsGood(j) = find(TempData(:,j),1,'last');
        TempData(IdsGood(j)+1:end,j) = NaN;
    end
    AllData.EMG_M{i,1}(:,1) = mean(TempData,2,'omitnan'); % 平均反应
    AllData.EMG_M{i,1}(:,2) = std(TempData,0,2,'omitnan')/sqrt(size(TempData,2)); % SEM
    AllData.EMG_Time{i,1} = AllEMGraw{1, 1}.interval*(1:size(TempData,1))'; % 时间（s）
    AllData.EMG{i,ROI_No} = TempData;
    
    TempData = AllData.Run{i,ROI_No};
    TrialNum = size(TempData,2);
    % 将末尾的连续0置为NaN
    for j = 1:TrialNum
        IdsGood(j) = find(TempData(:,j),1,'last');
        TempData(IdsGood(j)+1:end,j) = NaN;
    end
    AllData.Run_M{i,1}(:,1) = mean(TempData,2,'omitnan'); % 平均反应
    AllData.Run_M{i,1}(:,2) = std(TempData,0,2,'omitnan')/sqrt(size(TempData,2)); % SEM
    AllData.Run_Time{i,1} = AllRunraw{1, 1}.interval*(1:size(TempData,1))'; % 时间（s）
    AllData.Run{i,ROI_No} = TempData;
end
FileName = fullfile(SaveFolder,[TrialTitle,'AllData.mat']);
save(FileName,'AllData','RspOfSingleTrial');

disp('All Finished!');