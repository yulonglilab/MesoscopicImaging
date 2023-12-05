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