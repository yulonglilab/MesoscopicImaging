%% Spike2MultiInfoExtractVer7a
% Running time:20230629
%   ======================================================================================
% Fei Deng,20201210,增加了对于最后发光一段记录图片数目是否完整的检查，不完整则不进行提取
% Fei Deng,20201209,可以用于两个激
% Fei Deng,20201117,用于提取spike2记录的数据
% Fei Deng,20211206,用于提取spike2记录的数据
% Fei Deng,20211214,用于提取行为摄像头和spike2记录的数据
% Fei Deng,20211222,完整运行确保无误
% Fei Deng,20220117,增加代码自动计算输出视频时长和保存成像时刻信息
% Fei Deng,20220228,修正Spike2MultiInfoExtractVer6中保存response结果时错误替换的问题
% Fei Deng,20220313,Ratio改为Cor
% Spike2MultiInfoExtractVer6d,Fei
% Deng,20220510,修改视频图片读取方式，仅读取待显示帧,适用于总帧数很多，间隔较大的情形
% Fei Deng,20230601,Ver7，增加了输出荧光成像的视频(并根据不同状态添加文字标注），并直接指定视频的加速倍数
% Fei Deng,20230601,Ver7a，修复了小bug
% =======================================================================================
%   数据说明：
%   时间单位s
% =======================================================================================
close all;
clearvars -except CropList,clc;
% clear,clc,

ImCropListIJ = [23, 44, 311, 245]; % 成像视频截取区域，复制自ImageJ的坐标
ImCropList = ImCropListIJ+1;
ImCropList(3) = ImCropList(2);
ImCropList(2) = ImCropList(1)+ImCropListIJ(3);
ImCropList(4) =  ImCropListIJ(2)+ImCropListIJ(4)+1;

% 时间段截取
% 1.直接输入截取图片的时间范围
% xlimTime = [234.4 290.8]; % unit xTimeUnitList{xTimeUnitIdx}
% xlimTime = [234.4+3016/60 234.4+3235/60]; % unit xTimeUnitList{xTimeUnitIdx} 1340-1820
% xTimeUnitList = {'s','min','h'};
% xTimeUnitIdx = 2; % 第几种单位
% xlimTimeInSec = xlimTime*60^(xTimeUnitIdx-1); % unit xlimTimeInSec: second

% 2.直接输入截取图片的帧范围
BlockListS = [21641,23201,28801];
BlockListE = [22100,23800,29200];
BlockFrameNumList = cat(1,BlockListS,BlockListE)'; % 待截取block的帧数范围
% BlockFrameNumList =[9821,10080;17161,17380;17601,17900;18101,18330;18401,18710;18831,19100;]; % 待截取block的帧数范围
% BlockFrameNumS = 1251; % 待截取block的帧数范围
% BlockFrameNumE = 1450; % 待截取block的帧数范围
pseudoscaleList = [-0.2,2.5;
    -0.2,3]; % 成像的伪彩图显示scale

VideoFlag = 1; % 是否截取视频
VideoImFlag = 1;  % 是否截取成像视频
LabelFlag = 1; % 是否添加文字标注
LabelFlagStatus = 0; % 是否添加此时成像图片所对于状态的文字标注
TimeFlag = 1; % 视频标注中时间显示模式：1- s, 2-MM:SS, 3-HH:MM:SS
FrameInterNum = 1; % 生成response video的取图片间隔，每FrameInterNum取一帧

FrameIntvlStep = 5; % 输出视频时间尺度上间隔多少帧取一张，用于减小输出文件尺寸
ImFrameIntvlStep = 1; % 成像的输出视频时间尺度上间隔多少帧取一张，用于减小输出文件尺寸
BrightIndex = [0.8,1]; % 视频图像亮度调整系数
SpeedUpIdx = 20; % 输出的response video的加速倍数

TargetLabel = 'RspOfSession.mat'; % 输入目标文件夹/文件特征字符串
StatusLabels = {'REM','Wake','NREM'}; % Analysis by AccuSleep (1 = REM sleep, 2 = wakefulness, 3 = NREM sleep, 4 = undefined)
StatusColors = {[255,128,255]/255,[130,208,255]/255,[240,240,240]/255};
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
[file,path] = uigetfile('*RspOfSession.mat','请选择 RspOfSession.mat文件'); %文件路径
ChTagIdx = strfind(file,TargetLabel)-2;
savefolder = [path,'\Synchronized data'];
mkdir(savefolder);
TargetFileName = file;
TargetFileName(ChTagIdx) = '*';
cd(path);
targetPack = dir(TargetFileName); %待整理文件列表获取
ChNum = size(targetPack,1); % 待转换的帧数编号，自动计算
FileName = fullfile(path,file);
load(FileName);
ImagingTimeM = RspOfSession.TimeInSec;
ImPeriodT = (ImagingTimeM(end)-ImagingTimeM(1))/(length(ImagingTimeM)-1);

if exist('xlimTimeInSec')
    BlockListS = ceil(xlimTimeInSec(1)/ImPeriodT);
    BlockListE = floor(xlimTimeInSec(end)/ImPeriodT);
end
BlockFrameNumList = cat(1,BlockListS,BlockListE)'; % 待截取block的帧数范围

CropList = [];
%% 提取每个block的数据
for fi = 1:size(targetPack,1)
    cd(path);
    FileName = targetPack(fi).name;
    ChTag = FileName(ChTagIdx);
    FilePathName = fullfile(path,FileName);
    load(FilePathName);
    
    switch ChTag
        case 'G'
            TargetLabel = '_BVcor_Rsp.tif'; % 输入目标文件夹/文件特征字符串
            ChTagNo = 1; % 相机通道编号
        case 'R'
            TargetLabel = '_GVcor_Rsp.tif'; % 输入目标文件夹/文件特征字符串
            ChTagNo = 2; % 相机通道编号
    end
    TargetFileName = [FileName(1:strfind(FileName,'OfSession.mat')-5),TargetLabel];
    if LabelFlagStatus
        Idx = strfind(TargetFileName,'_Rsp.tif');
        FileImStatusG = [TargetFileName(1:Idx-1),'_mean images of status.mat'];
        load(FileImStatusG);
    end
    
    for bi = 1:size(BlockFrameNumList,1)
        BlockFrameNumS = BlockFrameNumList(bi,1); % 待截取block的起始帧数
        BlockFrameNumE = BlockFrameNumList(bi,2); % 待截取block的终止帧数
        ImagingTimeMb = ImagingTimeM(BlockFrameNumS:BlockFrameNumE)-ImagingTimeM(BlockFrameNumS);
        %     TableSync = {'Time (s)','Time (min)','Time (h)'};
        %     ImagingTimeMb_hms = repmat(ImagingTimeMb,1,3)./[1,60,3600];
        %         TableSync(2:size(ImagingTimeMb,1)+1,:) = num2cell(ImagingTimeMb_hms);
        %     filename = fullfile(savefolder,['ImagingTime_im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'.csv']);
        %     writecell(TableSync,filename);
        RspOfSessionBlock = cat(2,RspOfSession.Cor{1, 1}(BlockFrameNumS:BlockFrameNumE),RspOfSession.ChS{1, 1}(BlockFrameNumS:BlockFrameNumE),RspOfSession.ChNorm{1, 1}(BlockFrameNumS:BlockFrameNumE));
        TableSync = {'Time (s)','Cor','ChS','ChNorm'};
        RspOfSessionBlock = cat(2,ImagingTimeMb,RspOfSessionBlock);
        TableSync(2:size(ImagingTimeMb,1)+1,:) = num2cell(RspOfSessionBlock);
        filename = fullfile(savefolder,[FileName(1:strfind(FileName,'OfSession.mat')-1),'_',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'.csv']);
        disp(['Saving file: ',filename]);
        writecell(TableSync,filename);
        disp('File saved!');
        %% 截取成像视频
        if VideoImFlag == 1  % 是否截取视频
            disp(['Importing image series: ',TargetFileName]);
            % 读取图片
            fullFileName = fullfile(targetPack(fi).folder,TargetFileName);
            FileInfo = imfinfo(fullFileName);
            ImSz = [FileInfo(1).Height,FileInfo(1).Width,size(FileInfo,1)]; % 原始图片像素尺寸
            %     [ImSz1,ImSz852,ImSz3] = [FileInfo(1).Width,FileInfo(1).Height,size(FileInfo,1)]; % 原始图片像素尺寸
            ImSz1 = FileInfo(1).Height;
            ImSz2 = FileInfo(1).Width;
            FrameStart = BlockFrameNumS;
            FrameEnd = BlockFrameNumE;
            ImSz3 = FrameEnd-FrameStart+1;
            ImSz3 = ceil(ImSz3/ImFrameIntvlStep);
            IMstack = single(zeros(ImSz1,ImSz2,ImSz3));
            j = 0;
            for Idx = FrameStart:ImFrameIntvlStep:FrameEnd
                j = j+1;
                if mod(j,100)==0
                    disp(j);
                end
                IMstack(:,:,j) = imread(fullFileName,Idx);
            end
            
            if exist('ImCropList')
                ImCropX = ImCropList(1:2); % 读取已经暂存的每个行为学相机的ROI
                ImCropY = ImCropList(3:4);
                framesCrop = IMstack(ImCropY(1):ImCropY(2),ImCropX(1):ImCropX(2),:);
                [Sz1 Sz2 Sz3] = size(framesCrop);
            else
                framesCrop = IMstack;
            end
            
            % sigLims = [min(abs(wvcfs(:))), max(abs(wvcfs(:)))];
            sigLims = pseudoscaleList(ChTagNo,:);
            
            % Create video file
            fig = figure(105);
            figScale = 1;
            set(figure(fig),'Position',[10,50,Sz2*figScale,Sz1*figScale]);
            axis equal off;
            hold on;
            Idx = strfind(TargetFileName,'.tif');
            vidName = [TargetFileName(1:Idx-1),'_Frame',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE)];
            vidTitle = strcat(vidName,'_',datestr(now,'yymmddHHMM'),'_Min',num2str(sigLims(1)),'Max',num2str(sigLims(2)),'_',num2str(SpeedUpIdx),'X.avi');
%             vidObj = VideoWriter(vidTitle);
            vidObj = VideoWriter(fullfile(targetPack(fi).folder,vidTitle));
            vidFps = 1/ImPeriodT/ImFrameIntvlStep*SpeedUpIdx; % 视频帧率
            vidObj.FrameRate = vidFps;
            VideoT = ImSz3/vidObj.FrameRate;
            disp('Video time (s): ');
            disp(VideoT);
            
            open(vidObj);
            
            posX1 = round(0.3*Sz2);
            posY1 = round(0.055*Sz1);
            posX2 = round(0.95*Sz2);
            posY2 = round(0.055*Sz1);
            
            for j = 1:ImSz3
                %                 disp(j);
                clf;
                imshow(framesCrop(:,:,j),sigLims,'border','tight','initialmagnification','fit');
                colormap(gca, 'jet');
                if LabelFlag
                    PeriodOfVideoW = ImFrameIntvlStep*ImPeriodT;
                    box_color = {'red','green','yellow'};
                    FrameCurTime = PeriodOfVideoW*j;
                    switch TimeFlag
                        case 1
                            text_str = [num2str(round(FrameCurTime)),' s'];
                        case 2
                            text_str = datestr(FrameCurTime/3600/24,'MM:SS');
                        case 3
                            text_str = datestr(FrameCurTime/3600/24,'HH:MM:SS');
                    end
                    text(posX1,posY1,text_str,'Color',[0.99,0.99,0.99],'FontSize',20,'HorizontalAlignment','right');
                    % 添加状态标注
                    if LabelFlagStatus
                        IMofStatusBlock = IMofStatus.StatusStamp(FrameStart:ImFrameIntvlStep:FrameEnd);
                        StatusT = StatusLabels{IMofStatusBlock(j)};
                        text(posX2,posY2,StatusT,'Color',StatusColors{IMofStatusBlock(j)},'FontSize',20,'HorizontalAlignment','right');
                    end
                end
                writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
            end
            % Close video file
            close(vidObj);
        end
        
        if fi == 1
            %         RspOfSession.StartTimeInSec = 1.69; % !!!!!!!!!!!!!!!! 手动添加
            BlockStartT = ImagingTimeM((BlockFrameNumS))+RspOfSession.StartTimeInSec; % 所截取部分在整个session的起始时刻（s）
            BlockEndT = ImagingTimeM((BlockFrameNumE))+RspOfSession.StartTimeInSec; % 所截取部分在整个session的终止时刻（s）
            EEGb = EEG;
            EEGb.values = EEG.values(max([ceil(BlockStartT/EEG.interval),1]):min([ceil(BlockEndT/EEG.interval),EEG.length]));
            EEGb.length = length(EEGb.values);
            
            EMGb = EMG;
            EMGb.values = EMG.values(max([ceil(BlockStartT/EMG.interval),1]):min([ceil(BlockEndT/EMG.interval),EMG.length]));
            EMGb.length = length(EMGb.values);
            
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
            
            Timeb = Run.interval*(1:Runb.length)';
            TableSync = {'Time (s)','Velocity (cm/s)'};
            % Timeb = Timeb/60;
            % TableSync = {'Time (min)','Velocity (cm/s)'};
            Rundata = cat(2,Timeb,Velocity);
            TableSync(2:Runb.length+1,:) = num2cell(Rundata);
            filename = fullfile(savefolder,['Velocity_im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'.csv']);
            disp(['Saving file: ',filename]);
            writecell(TableSync,filename);
            disp('File saved!');
            
            Timeb = EEGb.interval*(1:EEGb.length)';
            TableSync = {'Time (s)','EEG (uV)','EMG (uV)'};
            % Timeb = Timeb/60;
            % TableSync = {'Time (min)','EEG (μV)','EMG (μV)'};
            EEGEMGdata = cat(2,Timeb,EEGb.values,EMGb.values);
            TableSync(2:EEGb.length+1,:) = num2cell(EEGEMGdata);
            filename = fullfile(savefolder,['EEGEMG_im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'.csv']);
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
            filename = fullfile(savefolder,['EEGEMG_im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'_Bin Step=',num2str(BinStep),'.csv']);
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
                    %                     VideoT = VideoTimeList(bi);
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
                                'FontSize',20,'BoxColor',box_color,'BoxOpacity',0,'TextColor','white');
%                             framesCropNoted(:,:,:,j) = insertText(framesCrop(:,:,:,j)*BrightIndex(i),position,text_str,...
%                                 'FontSize',30,'BoxColor',box_color,'BoxOpacity',0,'TextColor','white');
                            %     figure,imshow(framesCropNoted(:,:,:,j));
                        end
                    end
                    %% 输出视频
                    VideoNameW = fullfile(savefolder,[VideoName(1:strfind(VideoName,'mp4')-2),'_im',num2str(BlockFrameNumS),'-',num2str(BlockFrameNumE),'_T',num2str(frameStartT),'-',num2str(frameEndT),'s-',num2str(FrameIntvlStep),'_',num2str(SpeedUpIdx),'X.avi']);
                    VideoW = VideoWriter(VideoNameW,'Motion JPEG AVI');
                    VideoW = VideoWriter(VideoNameW);
                    % VideoW = VideoWriter(VideoNameW,'Uncompressed AVI');
                    % VideoW.Quality  = 100;
                    VideoW.FrameRate = VideoR.FrameRate/FrameIntvlStep*SpeedUpIdx; % 视频帧率
                    %                     VideoW.FrameRate = round(Sz4/VideoT);
                    VideoT = Sz4/VideoW.FrameRate;
                    disp('Video time (s): ');
                    disp(VideoT);
                    
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
disp('All Finished!');