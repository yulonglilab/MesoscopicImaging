%% Spike2MultiInfoExtractVer5
% Running time:20211222
%   ======================================================================================
% Fei Deng,20211222,用于进行大视场成像帧数和spike2记录的时间对应转换
%   =======================================================================================
FrameNoList = [2301/3,2550/3;
               2751,3000;
               3501,4200]; % 待转换的帧数编号，需要手动输入
SeqChNum = 2; % 待转换的帧数编号，需要手动输入
TypeFlag = 1;
%%
TargetFolder = uigetdir('','请选择指定的StimRcd文件夹');
cd(TargetFolder);
filetype = 'mat';
% [file,path] = uigetfile(['*.',filetype]); %文件路径
file = dir(['*.',filetype]); %文件路径
FileName = file.name;
load(FileName);
%%
FrameTimeList = cell(size(FrameNoList));
for fri = 1:size(FrameNoList,1)
    for fci = 1:size(FrameNoList,2)
        FrameNo = FrameNoList(fri,fci);
        switch TypeFlag
            case 1
                TimeSeconds = mean(Imaging.times([(FrameNo-1)*SeqChNum*2+1,FrameNo*SeqChNum*2]));
                TimeH = floor(TimeSeconds/3600);
                TimeM = floor((TimeSeconds-TimeH*3600)/60);
                TimeS = TimeSeconds-TimeH*3600-TimeM*60;
                TimeInSpike2HMS = [num2str(TimeH),':',num2str(TimeM),':',num2str(TimeS)];
                FrameTimeList{fri,fci} = TimeInSpike2HMS;
                disp(['Frame: ',num2str(FrameNo),' => Time (h:m:s): ',TimeInSpike2HMS]);
        end
    end
end
disp('Finished!');
%%
 TimeInHMS = '0:26:52';
  TimeInHMS = '0:26:54';
 Idx = find(TimeInHMS==':');
TimeInSec = str2num(TimeInHMS(1:Idx(1)-1))*3600+str2num(TimeInHMS(Idx(1)+1:Idx(2)-1))*60+str2num(TimeInHMS(Idx(2)+1:end));
VideoR.NumFrames/VideoR.FrameRate
% TimeConvert_Ver0(TypeFlag,Input,SeqChNum)
%% =======================================================================================
%% 
% function TimeConvert_Ver0(TypeFlag,Input,SeqChNum)
% switch TypeFlag
%     case 1
%         TimeSeconds = mean(Imaging.times([(FrameNo-1)*SeqChNum+1,FrameNo*SeqChNum]));
%         TimeH = floor(TimeSeconds/3600);
%         TimeM = floor((TimeSeconds-TimeH*3600)/60);
%         TimeS = TimeSeconds-TimeH*3600-TimeM*60;
%         disp([num2str(TimeH),':',num2str(TimeM),':',num2str(TimeS)]);
% end
% end