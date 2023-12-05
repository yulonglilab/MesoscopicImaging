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
