%% plotEEGEMGpatchSleepStateRspDiff_Ver1b
% Running time:20210625
%   =======================================================================================
% Fei Deng,20210622,for sleep wake analysis
% plotEEGEMGpatchSleepStateRspDiff_Ver1a,Fei Deng,20210625,增加EEG伪彩和RMS plotEMG的方式
% plotEEGEMGpatchSleepStateRspDiff_Ver1b,Fei Deng,20210626,仅plot需要的时间段
%   =======================================================================================
close all;
clear,clc;
colorpool.wake = [0.509, 0.815, 1];
colorpool.REM = [170, 170, 170]/255;
% colorpool.wake = [0.509, 0.815, 1];
% colorpool.REM = [1, 0.4, 0.8];
% colorpool.Ch=[0.01,0.631,0.01]; % Trace color of ROIs
colorpool.Ch = {[0,147,0]/255;
    [255,0,0]/255}; % Trace color of ROIs
YlimResponse = [-0.05,0.4;-0.05,0.05];
xTimeUnitList = {'s','min','h'};
% xTimeUnitIdx = 2;
xTimeUnitIdx = 2;
xTimeUnit = xTimeUnitList{xTimeUnitIdx};
% xlimTime = [0.516 1.1]; % unit xTimeUnitList{xTimeUnitIdx}
% xlimTime = [3.9 4.84]*60; % unit xTimeUnitList{xTimeUnitIdx}
xlimTime = [234.4 290.8]; % unit xTimeUnitList{xTimeUnitIdx}
% xlimTime = [0.525*60 1.1*60]; % unit xTimeUnitList{xTimeUnitIdx}
xticksTimeIntvl = 5; % unit xTimeUnitList{xTimeUnitIdx}
% xticksTimeIntvl = 1; % unit xTimeUnitList{xTimeUnitIdx}
% Image_f = 1; % 图像采样频率
EEGEMG_f = 1000; % EEGEMG采样频率
% CameraNum = 2;
MfileDir = pwd;
addpath(genpath(MfileDir));
ParentFolder = fileparts(MfileDir);
%%
% clearvars -except newpath;
% load EEG EMG data
SleepStateFolder = uigetdir('','请选择指定的SleepState文件夹');
cd(SleepStateFolder);
load('EEG_Accuformat.mat');
load('EMG_Accuformat.mat');
load('SleepStateLabel.mat');

% load('EMG_Accuformat.mat');
filetype = 'mat';

% load response data
TargetLabel = 'RspOfSession.mat'; % 输入目标文件夹/文件特征字符串
[file,path] = uigetfile(['*',TargetLabel],'请选择 RspOfSession.mat文件'); %文件路径
SaveFolder = [path,'Response'];
mkdir(SaveFolder);

labelsLen = size(labels,1);
stage = labels(1:labelsLen)';
% EEG and EMG parameter
% SampleIntervalE = EEG0.interval;
SampleIntervalE = 1/EEGEMG_f;
LenEEG = size(EEG,1);
SampleTotaltime = LenEEG*SampleIntervalE;
xtimeEEGall = (0:SampleIntervalE:SampleTotaltime-SampleIntervalE)/60^(xTimeUnitIdx-1); % 时间单位 unit xTimeUnitList{xTimeUnitIdx}
% xtimeEEG = (0:SampleInterval:SampleTotaltime-SampleInterval); % 时间单位s
if ~exist('xlimTime')
    xlimTime = [xtimeEEGall(1) xtimeEEGall(end)]; % unit hplot(xtimeEEG,EEG,'k');
end
xticksTime = ([xlimTime(1):xticksTimeIntvl:xlimTime(2)]);
% 生成EEG伪彩
Idx1 = find(xtimeEEGall<=xlimTime(1),1,'last');
Idx2 = find(xtimeEEGall>=xlimTime(end),1,'first');
params.Fs = EEGEMG_f;
params.fpass = [0 20];
params.tapers=[2 3];
movingwin=[EpochTime 0.1];
% data_YB0 = EEG;
% [spect_YB0, stimes0, sfreqs0]=mtspecgramc(data_YB0, movingwin, params);
% data_YB = EEG(Idx1-EpochTime*params.Fs:Idx2+EpochTime*params.Fs);
data_YB = EEG(Idx1:Idx2);
[spect_YB, stimes, sfreqs]=mtspecgramc(data_YB, movingwin, params);
stimes = stimes/60^(xTimeUnitIdx-1)+xtimeEEGall(Idx1);   %  统一时间单位 unit xTimeUnitList{xTimeUnitIdx}
% Take a median of spectrogram across 10 subjects
med_spect_YB = median(spect_YB,3);
PowerMap = pow2db(med_spect_YB(:,:)');
figure,
% imagesc(stimes/60/60, sfreqs, pow2db(med_spect_YB(:,:)'));
imagesc(stimes, sfreqs, PowerMap);
% imagesc(xtimeEEG(1:10000), sfreqs, pow2db(med_spect_YB((1:10000),:)'));
set(gca,'clim',[-20 50]);% set colorbar scale
axis xy;
ylabel('Frequency (Hz)','fontsize',13);
xlabel(['Time (',xTimeUnit,')'],'fontsize',13);
colormap('jet');
c = colorbar;
ylabel(c,'Power (dB)');
ylim(params.fpass);
% xlim([min(xtimeEEG) max(xtimeEEG)]);
set(gca,'TickDir','out');


ROI_PlotNum = 1;%ROI_Num-1;
ChTagIdx = strfind(file,TargetLabel)-2;
TargetFileName = file;
TargetFileName(ChTagIdx) = '*';
cd(path);
targetPack = dir(TargetFileName); %待整理文件列表获取

Nh = 2+ROI_PlotNum;
Nw = 2; % 多出一列用于画colorbar
gap = [.02 .03];
marg_h = [.1 .01];
marg_w = [.055 .1];
SzRatio_h = [2,1,1];
SzRatio_w = [30,1];

SreenSz = get(0,'ScreenSize');
%%
CameraNum = size(targetPack,1);
% plot single channel
for fi = 1:CameraNum
    %%
    FileName = targetPack(fi).name;
    ChTag = FileName(ChTagIdx);
    FilePathName = fullfile(path,FileName);
    load(FilePathName);
    SessionTitle{fi,1} = FileName(1:strfind(FileName,TargetLabel)-2);
    
    % h = figure();				% 创建图形窗口
    % h.Name = 'Sleep-Wake status&Response'; % 重命名图形窗口
    % hold on;
    %     % 最大化图形窗口，参考自： 原文链接：https://blog.csdn.net/am290333566/article/details/84581313
    %     warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	% 关闭相关的警告提示（因为调用了非公开接口）
    %     jFrame = get(h,'JavaFrame');	% 获取底层 Java 结构相关句柄吧
    %     pause(0.1);					% 在 Win 10，Matlab 2017b 环境下不加停顿会报 Java 底层错误。各人根据需要可以进行实验验证
    %     set(jFrame,'Maximized',1);	%设置其最大化为真（0 为假）
    %     pause(0.1);					% 个人实践中发现如果不停顿，窗口可能来不及变化，所获取的窗口大小还是原来的尺寸。各人根据需要可以进行实验验证
    %     warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		% 打开相关警告设置
    figure('name','Sleep-Wake status&Response');
    set(gcf,'Position',[20 400 1500 300]);
%     set(gcf,'Position',[20 400 1500 500]);
    ha = tight_subplot_Ver1(Nh, Nw, gap, marg_h, marg_w,SzRatio_h,SzRatio_w);
    % ha = tight_subplot_Ver1(2+ROI_PlotNum,1,[.02 .03],[.1 .01],[.055 .01],[2,1,1],[1,1,1]);
    
    %%
    % plot response
    for i = 1:ROI_PlotNum
        %     IdxROI = i+1;
        SampleNum = size(RspOfSession.TimeInSec,1);
        %     SampleInterval = ImPeriod;
        %     SampleTotaltime = SampleNum*SampleInterval;
        switch xTimeUnitIdx
            case 1
                xtimeIM = RspOfSession.TimeInSec+RspOfSession.StartTimeInSec; % 时间单位s
            case 2
                xtimeIM = RspOfSession.TimeInMin+RspOfSession.StartTimeInSec/60; % 时间单位min
            case 3
                xtimeIM = RspOfSession.TimeInHrs+RspOfSession.StartTimeInSec/3600; % 时间单位h
        end
        % xtimeIM = (0:SampleInterval:SampleTotaltime-SampleInterval); % 时间单位s
        Idx1 = find(xtimeIM<=xlimTime(1),1,'last');
        Idx2 = find(xtimeIM>=xlimTime(end),1,'first');
        Idx_axes = (i-1)*Nw+1;
        axes(ha(Idx_axes));
        h_axes{Idx_axes} = plot(xtimeIM(Idx1:Idx2),RspOfSession.Cor{1, 1}(Idx1:Idx2),'color',colorpool.Ch{fi});
        h_axes{Idx_axes}.LineWidth = 1; 
        set(gca,'TickDir','out');
        % axis tight
        set(gca,'box','off');
        set(gca,'color','none');
        ha(i).XAxis.Visible = 'off';
        xlim(xlimTime);
        ylim(YlimResponse(fi,:));
        ylabel(' \DeltaF/F_{0}');
        %     ylabel([SingleROI.Title{IdxROI,1},' \DeltaF/F_{0}']);
        hold on;
        % label special points
%         ImIdx = [14115,14266,14396,14468,14792,15019,15370,15820,17173,17216,17384];
        ImIdx = [14115,14266,14468,14792,15463,15820,16008,16707,16945,17384];
          Idx_axes = (i-1)*Nw+1;
        axes(ha(Idx_axes));
        scatter(xtimeIM(ImIdx),RspOfSession.Cor{1, 1}(ImIdx),15,'MarkerEdgeColor','r');
    end
    
    % plot EEG EMG
    %     % EEG, trace
    %     axes(ha(ROI_PlotNum+1));
    %     plot(xtimeEEG,EEG,'k');
    %     set(gca,'TickDir','out');
    %     % axis tight
    %     set(gca,'box','off');
    %     set(gca,'color','none');
    %     ha(ROI_PlotNum+1).XAxis.Visible = 'off';
    %     xlim(xlimTime);
    %     ylim([-500 500]);
    %     ylabel('EEG (uV)');
    
    % EEG, 时频图
    Idx_axes = ROI_PlotNum*Nw+1;
    axes(ha(Idx_axes));
    h_axes{Idx_axes} = imagesc(stimes, sfreqs, PowerMap);
    set(gca,'clim',[-10 50]);% set colorbar scale
    %     set(gca,'clim',[min(PowerMap(:))-abs(min(PowerMap(:)))*0.5,max(PowerMap(:))*1.5]);% set colorbar scale
    axis xy;
    ylabel('Frequency (Hz)');
    colormap('jet');
    %     colorbar('eastoutside')
    %     colorbar('off')
    %     c = colorbar;
    %     ylabel(c,'Power (dB)');
    colormapEEG = colormap;
    ylim(params.fpass);
    xlim(xlimTime);
    set(gca,'TickDir','out');
    set(gca,'box','off');
    set(gca,'color','none');
    ha(ROI_PlotNum*Nw+1).XAxis.Visible = 'off';
    
    % EMG
    %         EMGf = highpass(EMG,30,EEGEMG_f);
    %         xtimeEMGf = xtimeEEG;
    EpochRMS = 1; % RMS计算时平均的窗口时长 (s)
    Idx1 = find(xtimeEEGall<=xlimTime(1),1,'last');
    Idx2 = find(xtimeEEGall>=xlimTime(end),1,'first');
    % xtimeEMG = xtimeEEGall(Idx1:Idx2);
    
    xtimeEMGf = xtimeEEGall(Idx1):EpochRMS/60^(xTimeUnitIdx-1):xtimeEEGall(Idx2);
    EMGf = xtimeEMGf;
    for i = 1:size(xtimeEMGf,2)
        if mod(i,100)==0
            disp(i)
        end
        IdxT = min(Idx1+EpochRMS*EEGEMG_f-1,Idx2);
        EMGf(i) = rms(EMG(Idx1:IdxT));
        Idx1 = IdxT+1;
    end
    %       EMGf2 = highpass(EMGf,30,EEGEMG_f);
    Idx_axes = (ROI_PlotNum+1)*Nw+1;
    axes(ha(Idx_axes));
    %     figure,
    h_axes{Idx_axes} = plot(xtimeEMGf,EMGf,'k');
     h_axes{Idx_axes}.LineWidth = 1; 
    set(gca,'TickDir','out');
    % axis tight
    set(gca,'box','off');
    set(gca,'color','none');
    % ha(ROI_PlotNum+2).XAxis.Visible = 'off';
    xlim(xlimTime);
    %     ylim([-1000 1000]);
    %     ylabel('EMG (uV)');
    ylim([prctile(EMGf,1)*0.9 prctile(EMGf,99.95)]);
    %     ylim([34 135]);
    ylabel('EMG RMS (uV)');
    xticks(xticksTime);
    xlabel(['Time (',xTimeUnit,')']);
    %     delete(h_axes{Idx_axes});
    %% Wake patch （根据AccuSleep编号，2表示wake）
    stage2 = [0,stage];
    TargetStage = zeros(size(stage2));
    TargetStage(find(stage2==2)) = 1; %（根据AccuSleep编号，2表示wake）
    TargetStageDiff = diff(TargetStage);
    TargetStageOn = xtimeEEGall((find(TargetStageDiff==1)-1)*EpochTime/SampleIntervalE+1);
    TargetStageOff = xtimeEEGall((find(TargetStageDiff==-1)-1)*EpochTime/SampleIntervalE);
    TrialCount = length(TargetStageOn);
    ticky = get(gca,'ylim');
    hold on;
    Idx1 = find(TargetStageOn<=xlimTime(1),1,'last');
    Idx2 = find(TargetStageOn>=xlimTime(end),1,'first');
    %     for i = 1:TrialCount
    k = 0;
    for i = Idx1:Idx2
        k = k+1;
        h1(k) = area([TargetStageOn(i) TargetStageOff(i)], [ticky(2) ticky(2)],ticky(1), 'FaceColor',colorpool.wake,'edgecolor', 'none');
        uistack(h1(k),'down');
    end
    
    %% REM patch （根据AccuSleep编号，1表示REM）
    stage2 = [0,stage,0];
    TargetStage = zeros(size(stage2));
    TargetStage(find(stage2==1))=1; %（根据AccuSleep编号，1表示REM）
    TargetStageDiff=diff(TargetStage);
    TargetStageOn=xtimeEEGall((find(TargetStageDiff==1)-1)*EpochTime/SampleIntervalE+1);
    TargetStageOff=xtimeEEGall((find(TargetStageDiff==-1)-1)*EpochTime/SampleIntervalE);
    TrialCount=length(TargetStageOn);
    %     ticky=get(gca,'ylim');
    hold on;
    Idx1 = find(TargetStageOn<=xlimTime(1),1,'last');
    Idx2 = find(TargetStageOn>=xlimTime(end),1,'first');
    %     for i = 1:TrialCount
    k = 0;
    for i = Idx1:Idx2
        k = k+1;
        h2(k) = area([TargetStageOn(i) TargetStageOff(i)], [ticky(2) ticky(2)],ticky(1), 'FaceColor',colorpool.REM,'edgecolor', 'none');
        uistack(h2(k),'down');
    end
    
    shading interp
    alpha(0.3)
    
    for i = Nw:Nw:Nw*Nh
        axes(ha(i));
        set(gca,'box','off');
        set(gca,'color','none');
        ha(i).XAxis.Visible = 'off';
        ha(i).YAxis.Visible = 'off';
    end
    axes(ha((ROI_PlotNum+1)*Nw));
    set(gca,'clim',[-20 50]);% set colorbar scale
%     colormap(colormapEEG); 
    colormap(CustomColormapEEGshow2); 
    colobarj = colorbar;
    ylabel(colobarj,'Power (dB)');
    colobarj.TickDirection='out';
    colobarj.Box='off';
    colobarj.Position(3) = 2*colobarj.Position(3);  % 改变系数0.3（设置合适的宽度）
    
    set(gcf,'Renderer', 'Painter');
%     FileName = fullfile(SaveFolder,[SessionTitle{fi,1},'_Status&Traces_ROI ',num2str(ROI_PlotNum),'_',num2str(xlimTime(1)),'-',num2str(xlimTime(2)),xTimeUnit,'.pdf']);
%     saveas(gcf,FileName);
%     FileName = fullfile(SaveFolder,[SessionTitle{fi,1},'_Status&Traces_ROI ',num2str(ROI_PlotNum),'_',num2str(xlimTime(1)),'-',num2str(xlimTime(2)),xTimeUnit,'.fig']);
%     saveas(gcf,FileName);
FileName = fullfile(SaveFolder,[SessionTitle{fi,1},'_Status&Traces_ROI ',num2str(ROI_PlotNum),'_',num2str(xlimTime(1)),'-',num2str(xlimTime(2)),xTimeUnit,'_label.pdf']);
saveas(gcf,FileName);
FileName = fullfile(SaveFolder,[SessionTitle{fi,1},'_Status&Traces_ROI ',num2str(ROI_PlotNum),'_',num2str(xlimTime(1)),'-',num2str(xlimTime(2)),xTimeUnit,'_label.fig']);
saveas(gcf,FileName);

    %% 计算均值
    % StatusStamp = imresize(uint8(stage') ,[SampleNum,1]); % 每帧图像对应的状态
    % StatusLabels = {'REM','Wake','NREM'};
    % StatusMeanIM = zeros(ImSz1,ImSz2,3,'single');
    % for i = 1:3
    %     % Calculate average response of each ROI
    %     IdxStatus = find(StatusStamp==i);
    %     StatusData = SingleROI.Rsp(:,IdxStatus);
    %     StatusMean(:,i) = mean(StatusData,2);
    %     StatusSEM(:,i) = std(StatusData,0,2)/sqrt(size(StatusData,2));
    %     % Generate average image of each status
    %     StatusMeanIM(:,:,i) = mean(IMseries(:,:,IdxStatus),3);
    %     IMsAve = StatusMeanIM(:,:,i);
    %     %         figure,imshow(IMsAve,[]);
    %     objt = Tiff([SaveFolder,'\',SessionTitle{fi,1},'_',StatusLabels{i},'.tif'],'w');
    %     tagstruct.ImageLength = size(IMsAve,1);
    %     tagstruct.ImageWidth = size(IMsAve,2);
    %     tagstruct.Photometric = Tiff.Photometric.MinIsBlack; %info.PhotometricInterpretation;  %
    %     tagstruct.BitsPerSample = info.BitDepth; %32;
    %     tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; %info.SampleFormat; %
    %     tagstruct.SamplesPerPixel = info.SamplesPerPixel; %1;
    %     tagstruct.RowsPerStrip    = info.RowsPerStrip;
    %     tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    %     tagstruct.Software = 'MATLAB';
    %     objt.setTag(tagstruct);
    %     objt.write(IMsAve);
    %     objt.close();
    % end
    % FileName = fullfile(SaveFolder,'Quantified Rsp of sleep states.mat');
    % save(FileName,'StatusMean','StatusSEM');
    
end
%% =======================================================================================
%% 函数tight_subplot，用于绘图
function [ha, pos] = tight_subplot_Ver1(Nh, Nw, gap, marg_h, marg_w,SzRatio_h,SzRatio_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end
if nargin<6 || isempty(SzRatio_h); SzRatio_h = ones(Nh); end
if nargin<7; SzRatio_w = ones(Nw); end

if numel(gap)==1;
    gap = [gap gap];
end
if numel(marg_w)==1;
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1;
    marg_h = [marg_h marg_h];
end

axh_unit = (1-sum(marg_h)-(Nh-1)*gap(1))/sum(SzRatio_h);
axw_unit = (1-sum(marg_w)-(Nw-1)*gap(2))/sum(SzRatio_w);

% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    axh = axh_unit*SzRatio_h(ih);
    if ih == 1
        py = 1-marg_h(2)-axh;
    else
        py = py-gap(1)-axh;
    end
    for iw = 1:Nw
        axw = axw_unit*SzRatio_w(iw);
        if iw == 1
            px = marg_w(1);
        else
            px = px+axw+gap(2);
        end
        
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
    end
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end