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