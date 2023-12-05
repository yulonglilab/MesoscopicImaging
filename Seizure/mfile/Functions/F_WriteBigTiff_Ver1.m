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