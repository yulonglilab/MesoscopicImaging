%% StimRspSingleTrial_Ver2a
% Running time:20220208
%   =======================================================================================
% Fei Deng,20210621,
%   =======================================================================================
%%
close all,
clear,clc;
%% Parameter setting
CameraNum = 2; % 使用的相机数目
StimTagList = {{'10Hz 10s','20Hz 10s','30Hz 10s','40Hz 10s','50Hz 10s','60Hz 10s','70Hz 10s','80Hz 10s','90Hz 10s','100Hz 10s'};
    {'50Hz 1s','50Hz 2s','50Hz 5s','50Hz 10s','50Hz 20s','50Hz 30s'};
    {'5Hz 10s','10Hz 10s','20Hz 10s','30Hz 10s','40Hz 10s','50Hz 10s'}}; % 刺激名称，每个session一行pseudoscaleMin = -0.05; % 伪彩图显示最小值,for G channel
pseudoscaleMax = 0.2; % 伪彩图显示最大值,for G channel
% pseudoscaleMin = -0.03; % 伪彩图显示最小值,for R channel
% pseudoscaleMax = 0.05; % 伪彩图显示最大值,for R channel
IMscale = 0.7; % Resize图片系数  IMscale = 0.7;
ImBinning = 4; % 实际成像时的binning数值  ImBinning = 3;
pixelSzBin1 = 2211.358; % 成像时binning为1X1时的标尺，pixels / cm， 20211031测定校准标尺
pixelsize = pixelSzBin1/ImBinning*IMscale;
compression = 0;
Image_T = 0.2; % 每帧图像周期（s）
% BaseDurList = [3;3]; %每个session的每个trial baseline计算时使用时间(s)
BaseDurList = ones(6,1)*10; %每个session的每个trial baseline计算时使用时间(s)

%%
MfileDir = pwd;
addpath(genpath(MfileDir));
ParentFolder = fileparts(MfileDir);
mainPath = uigetdir('','请选择指定的含TrialRsp.mat的文件夹');
cd(mainPath);

SaveFolder = [mainPath,'\RspOfSingleTrial'];

% mainPath = [ParentFolder,'\rawdata\Processing results'];
% mainPath = [ParentFolder,'\Processing results'];
% SaveFolder = [ParentFolder,'\Processing results\RspOfSingleTrial'];
mkdir(SaveFolder);
cd(mainPath);
TargetList = dir('*TrialRsp.mat');
TargetNum = size(TargetList,1);
SessionNum = TargetNum/CameraNum;


%% Load WholeCortexMask and ROIs in ImageJ
% load WholeCortexMask
MaskName = [SaveFolder,'\WholeCortexMask_reScale.tif'];
WholeCortexMask = single(imread(MaskName));
% WholeCortexMask0(find(WholeCortexMask0==0)) = nan;
% WholeCortexMask0(find(WholeCortexMask0==255)) = 1;
% WholeCortexMask = uint8(imresize(WholeCortexMask0,IMscale,'Method','nearest'));
%     figure,imshow(WholeCortexMask);
[ImSz1,ImSz2] = size(WholeCortexMask);
% FileName = [MaskName(1:end-4),'_reScale.tif'];% num2str(ImSz1),'pix',num2str(ImSz2),'pix.tif'
% obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
% WriteIMG(obj,WholeCortexMask'); % save averaged image
% close(obj);


%load ROIs
ROI_Name = [SaveFolder,'\RspRoiSet.zip'];
ROIs = ReadImageJROI(ROI_Name);
for ri = 1:size(ROIs,2)
    ROI = ROIs{ri};
    switch ROI.strType
        case 'Rectangle'
            ROIxy = ROI.vnRectBounds;
            mask = logical(zeros(ImSz1,ImSz2));
            mask(ROIxy(1):ROIxy(3),ROIxy(2):ROIxy(4)) = 1;
        case 'Polygon'
            mask = poly2mask(ROI.mnCoordinates(:, 1), ROI.mnCoordinates(:, 2),ImSz1,ImSz2);
        case 'Freehand'
            mask = poly2mask(ROI.mnCoordinates(:, 1), ROI.mnCoordinates(:, 2),ImSz1,ImSz2);
        otherwise
            disp('ROI type ERROR!');
    end
    [mask_r{ri}, mask_c{ri}] = find(mask > 0);
    mask_id{ri} = find(mask > 0);
    figure,
    imshowpair(mask,WholeCortexMask,'falsecolor'),title('After Registration');
    %     title(ROI.strName);
    title(['ROI ',num2str(ri)]);
end
close all;
%% Load response mat file and Analysis of response of each trial
for fi = 1:SessionNum
    for ci = 1:CameraNum
        % Load response mat file
        TargetPath = TargetList(ci+(fi-1)*CameraNum).name; % 为方便记录，所有原始数据按照session放置
        %     TargetPath = TargetList{fi};
        %     [filepath,SessionTitle{fi,1},ext] = fileparts(fileparts(TargetPath));
        %     if ~isempty(ext) % 避免文件夹名称含小数点导致的fileparts函数识别出错的情况
        %         SessionTitle{fi,1} = cat(2,SessionTitle{fi,1},ext);
        %     end
        %     TempTitle = SessionTitle{fi,1};
        %     Idx = find(TempTitle=='_',1,'last');
        %     SessionTitle{fi,1} = [TempTitle(1:Idx-3),TempTitle(Idx:end)];
        SessionTitle{fi,1} = TargetPath(1:strfind(TargetPath,'_TrialRsp.mat')-1);
        load(TargetPath);
        StimModeNum = size(IMseriesTrialRspRatio,2); % 刺激模式总数
        StimModeCount = zeros(1,StimModeNum);   % 每种刺激模式总数，即对应trial数
        
        % Analysis of response of each trial
        RspOfSingleTrial.Ratio = cell(StimModeNum,size(ROIs,2)); % 比值通道,每种刺激模式(行）和每个ROI（列）进行组合
        for si = 1:StimModeNum
            IMseriesRsp = IMseriesTrialRspRatio{si};
            [ImSz1,ImSz2,TrialFrames,StimModeCount(si)] = size(IMseriesRsp);
            % response of each ROIs in each single trials
            % calculate response
            for ti = 1:StimModeCount(si) % 该类型trial总数
                RspRatio = IMseriesRsp(:,:,:,ti);
                %         figure,imshow(mean(RspRatio,3,'omitnan'),[]);
                %         colormap jet;
                %         figure,imshow(max(RspRatio,[],3,'omitnan'),[-0.03 0.2]);
                %         colormap jet;
                RspRatio = reshape(RspRatio,[],TrialFrames);
                %             RspRatio = RspRatio';
                for ri = 1:size(ROIs,2)
                    %                 RspOfSingleTrial.Ratio{si,ri} = cat(2,RspOfSingleTrial.Ratio{si,ri},mean(RspRatio(:,mask_id{ri}),2,'omitnan')); % 比值的反应
                    RspOfSingleTrial.Ratio{si,ri} = cat(1,RspOfSingleTrial.Ratio{si,ri},mean(RspRatio(mask_id{ri},:),1,'omitnan')); % 比值的反应
                end
            end
            figure,
            for ri = 1:size(ROIs,2)
                clf;
                xRange = (0:60/Image_T)*Image_T-BaseDurList(fi);
                TempData = RspOfSingleTrial.Ratio{si,ri};
                TempData = TempData(:,1:size(xRange,2));
                %             imagesc((0:TrialFrames-1)*Image_T-BaseDurList(fi),1:StimModeCount(si),RspOfSingleTrial.Ratio{si,ri},[pseudoscaleMin,pseudoscaleMax]);
                
                imagesc(xRange,1:size(TempData,1),TempData,[pseudoscaleMin,pseudoscaleMax]);
                colormap jet;
                colorbarNote = colorbar;
                colorbarNote.Label.String = '\DeltaF/F_0';
                colorbarNote.Label.Rotation = 270;
                colorbarNote.Label.FontSize = 16;
                colorbarNote.Ticks = [pseudoscaleMin,pseudoscaleMax];
                set(gca,'XTick',-10:10:50,'fontsize',14);
                %         set(gca,'XTick',[0,round(TrialFrames/5):round(TrialFrames/5):TrialFrames]*Image_T);
                xlabel('Time(s)','fontsize',16);
                set(gca,'YTick',1:1:StimModeCount(si),'fontsize',14);
                ylabel('Trial No.','fontsize',16);
                %             title('Heat Map of Single Trial','fontsize',20);
                %             title(['Heat Map of Single Trial, ROI ',num2str(ri)],'fontsize',20);
                title(['ROI ',num2str(ri)],'fontsize',20);
                saveas(gcf,[SaveFolder,'\',SessionTitle{fi,1},'_',StimTag{si},'_ROI',num2str(ri),'.fig']);
                saveas(gcf,[SaveFolder,'\',SessionTitle{fi,1},'_',StimTag{si},'_ROI',num2str(ri),'.jpg']);
            end
            %每个ROI的多个trial平均值
            for ri = 1:size(ROIs,2)
                Value = RspOfSingleTrial.Ratio{si,ri}; % 反应值
                Value = Value';
                ValueM = mean(Value,2); % 平均反应
                ValueSEM = std(Value,0,2,'omitnan')/sqrt(size(Value,2)); % SEM
                RspOfSingleTrial.Ratio_M{si,ri} = cat(2,ValueM,ValueSEM);
                RspOfSingleTrial.Ratio_Stim{si,1}(:,2*ri-1:2*ri) = RspOfSingleTrial.Ratio_M{si,ri}; % 按照刺激类型归类
            end
        end
        FileName = [SaveFolder,'\',SessionTitle{fi,1},'_RspOfSingleTrial.mat']
        save(FileName,'RspOfSingleTrial','-v7.3');
        clear RspOfSingleTrial;
        close all;
    end
end
disp('Finished!');