%% F_StatusMeanIM_Ver1
% Running time:20220407
%   ======================================================================================
% Fei Deng,20220407,用于划分SleepWakeState
% =======================================================================================
%   数据说明：
%   时间单位s
%   F_StatusMeanIM_Ver1，增加每种类型图片帧数和类型名称的记录
% =======================================================================================
function IMofStatus = F_StatusMeanIM_Ver1(IMseries3D,StatusStamp,StatusLabels,FileNameSuffix,pixelsize,compression)
%% 计算均值
[ImSz1 ImSz2 ImSz3] = size(IMseries3D);
StatusNum = length(unique(StatusStamp));
StatusMeanIM = zeros(ImSz1,ImSz2,3,'single');
IMofStatus.Label = StatusLabels;
IMofStatus.StatusStamp = StatusStamp;
for i = 1:StatusNum
    IdxStatus = find(StatusStamp==i);
    IMofStatus.FrameNumber(i) = length(IdxStatus);
    IMsAve = mean(IMseries3D(:,:,IdxStatus),3);
    IMofStatus.Mean(:,:,i)= IMsAve;
    IMofStatus.SEM(:,:,i) = std(IMseries3D(:,:,IdxStatus),0,3,'omitnan')/sqrt(length(IdxStatus)); % SEM
    FileName = [FileNameSuffix,StatusLabels{i},'.tif'];
    obj = Fast_BigTiff_Write(FileName,pixelsize,compression);
    WriteIMG(obj,IMsAve'); % save averaged image
    close(obj);
end
FileName = [FileNameSuffix,'mean images of status.mat'];
save(FileName,'IMofStatus');