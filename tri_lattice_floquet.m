%构造tri-lattice-array floquet
clear, clc
% R0=110;%建模的中心半径
% R1=40;%建模的波动半径
d=10;
%%%%%%%%%%%%a矩阵为输出旋转相位矩阵
units = 'mm';
Units = 'mm';
Type = 'csv';

% add paths to the required m-files.
relPath = 'E:\老电脑桌面2020\E文档\资料\matlab与hfss连接\MatlaHFSSApi/';
addpath(relPath);
hfssIncludePaths(relPath);
%%%%%%%bianliang

%%%%%%%%%脚本开始准备
false = 0;
true = 1;

fileName = ['hexagonal_element'];%文件名
% PrjFile = ['E:\资料\matlab与hfss连接\my_tmptext\',fileName,'.aedt'];
% tmpDataFile = 'E:\资料\matlab与hfss连接\my_tmptext\tmpData.m'
tmpScriptFile = ['E:\老电脑桌面2020\any OAM mode\shape_based script\',fileName,'.vbs'];%stript文件位置
% HFSS Executable Path.
% hfssExePath = '"D:\SimulationWare\HFSS19\AnsysEM19.0\Win64\ansysedt.exe"'; 
    fid = fopen(tmpScriptFile, 'wt');   % 'wt'表示以文本模式打开文件，可写，覆盖原有内容 
   % 创建一个新的工程并插入一个新的设计
    hfssNewProject(fid);
    hfssInsertDesign(fid, fileName);
%     hfssPre(fid);
    %%%%%%%%%%%%%design Variables
%     f=7;%GHz
%     d=10;%周期
    h=3;
    b1=2.6;
    b2=7.8;
    s1=0.3;
    s2=2.6;
    da=20;
%     g1=d-b1;
%     g2=d-b2;
    ob_h=100;%观测高度
%     hfssaddVar(fid, 'f0', f, []);
    hfssaddVar(fid, 'd', d, units);
    hfssaddVar(fid, 'h', h, units);
    hfssaddVar(fid, 'b1', b1, units);
    hfssaddVar(fid, 'b2', b2, units);  
    hfssaddVar(fid, 's1', s1, units);
    hfssaddVar(fid, 's2', s2, units);  
    hfssaddVar(fid, 'da', da, units);  
    hfssaddVar(fid, 'p', 'sqrt(3)*d', units);
    hfssaddVar(fid, 'angl', 0 , 'deg');
%     hfssaddVar(fid, 'N', N, []); 
%     hfssaddVar(fid, 'ob_h', ob_h, units); 
%     hfssaddVar(fid, 'g1', g1, units);
%     hfssaddVar(fid, 'g2', g2, units); 
    %%%%%%%%Array of metasrface(square loop)
%     squareloopName = cell(arraynum, 1);%方片名字共arraynum个
%     squareloopName1 = cell(arraynum, 1);
%     hexagonalName1 = cell(arraynum, 1);
%     numb=1;
%%  %%%%开始计算建模
    %空气盒子
    hfssBox(fid, 'airbox',{'-d/2','-p/2','-da'},{'d','p','h+2*da'}, Units);
    %设置sub和材料
    hfssBox(fid, 'sub',{'-d/2','-p/2','0'},{'d','p','h'}, Units);
    Er=2.65;
    bSigma=0;
    tanDelta=0.002;
    hfssAddMaterial(fid, 'myf4b', Er, bSigma, tanDelta)
    hfssAssignMaterial(fid, 'sub', 'myf4b')
    %金属地
    Start_gud = {'-d/2','-p/2','0'};
    hfssRectangle(fid, 'ground', 'Z', Start_gud , 'd', 'p', units);  
    CellGroundName={'ground'};
    hfssAssignPE(fid, 'gud', CellGroundName,  false);%PEC 

%% 建立单元patch啦
patch = cell(4,1);%原子的四个金属片
patch1 = cell(4,1);
for i=1:4
        %add local CS1
            CSName = ['CS',num2str(i)];
            if i==1
                  Origin = {'-d/2','0','h'};
            end
            if i==2
                  Origin = {'d/2','0','h'}; 
            end
            if i==3
                  Origin = {'0','p/2','h'};
            end
            if i==4
                  Origin = {'0','-p/2','h'};
            end
            XAxisVec = [1 0 0];
            YAxisVec = [0 1 0];
            hfssCreateRelativeCS(fid, CSName, Origin,XAxisVec, YAxisVec, units);
            hfssSetWCS(fid, CSName);%%跳到建立的相对坐标建模loop
  %%创建具体
  patch{i}=['patch',num2str(i)];
  patch1{i}=['patch1',num2str(i)];
   Start_squareloop = {0,0,0};
%   hfssEllipse(fid, patch{i}, 'Z', [0,0,0], 'b1', 'b2/b1', units) 
%   hfssEllipse(fid, patch1{i}, 'Z', [0,0,0], 's1', 's2/s1', units)
            Start_squareloop = {'-b1/2','-b2/2','h-h'};
            hfssRectangle(fid,  patch{i}, 'Z', Start_squareloop, 'b1', 'b2', units);
            Start_squareloop1 = {'-b1/2+s1','-b2/2+s2','h-h'};
            hfssRectangle(fid, patch1{i}, 'Z', Start_squareloop1, 'b1-2*s1', 'b2-2*s2', units);
  hfssSubtract(fid, patch{i}, patch1{i});
  CellR_patchName={patch{i}};
  hfssRotate_y(fid, CellR_patchName, 'Z', 'angl');
   hfssSetWCS(fid, 'Global');
   %%%%%结束创建六边形
end
%和到一块   
  hfssUnite2(fid, 4,'patch')

  %% cut掉 不需要的部分
  %1.建立被cut的空心盒子
hfssBox(fid, 'cut_airbox',{'-d/2','-p/2','-da'},{'d','p','h+2*da'}, Units);%仿真区域
hfssBox(fid, 'cut',{'-d/2-da','-p/2-da','-da'},{'d+2*da','p+2*da','h+2*da'}, Units);%将要被cut的空心盒子
CutName={'cut'};
CutairboxName={'cut_airbox'};
hfssSubtract(fid, CutName, CutairboxName);
%2.金属剪掉空心盒子
hfssSubtract(fid, patch{1}, 'cut');

%%      %%%%%%%边界条件
  CellpatchName={'patch1'};
hfssAssignPE(fid, 'patch',  CellpatchName,  false);%PEC


%% 设置主从边界
hfssAssignMaster_y(fid, 'Master1', '10', [-5, 8.66025403784439, -20], ...
	                 [-5, 8.66025403784439, 23], 'mm', false);
hfssAssignSlave_y(fid, 'Slave1', '12', [5, 8.66025403784439, -20], ...
	                 [5, 8.66025403784439, 23], 'mm','Master1', true);
hfssAssignMaster_y(fid, 'Master2', '9', [5, -8.66025403784439, -20], ...
	                 [5, -8.66025403784439, 23], 'mm', true);
hfssAssignSlave_y(fid, 'Slave2', '11', [5, 8.66025403784439, -20], ...
	                 [5, 8.66025403784439, 23], 'mm','Master2', false);

% hfssAssignPE(fid, 'hexagonalground', hexagonalName1,  false);
% hfssAssignPE(fid, 'PEsquareloop', squareloopName,  false);
%% 求解设置 
hfssInsertSolution(fid, 'Setup1', 9, 0.05, 20) %1（fid）,2(name),3(fGHz),4(maxDeltaS)可无，5（maxPass）可无
hfssInterpolatingSweep_y(fid,'Sweep','Setup1', 6, 21, 0.5);
%% 添加变量
% fprintf(fid, 'Set oModule = oDesign.GetModule("OutputVariable")\n');
% %PCR
% fprintf(fid, 'oModule.CreateOutputVariable "PCR",  _\n');
% fprintf(fid, ' "(mag(S(FloquetPort1:2,FloquetPort1:2)-S(FloquetPort1:1,FloquetPort1:1)))^2/((m" & _\n');
% fprintf(fid, ' "ag(S(FloquetPort1:2,FloquetPort1:2)-S(FloquetPort1:1,FloquetPort1:1)))^2+(mag(" & _\n');
% fprintf(fid, ' "S(FloquetPort1:2,FloquetPort1:2)+S(FloquetPort1:1,FloquetPort1:1)))^2)",  _\n');
% fprintf(fid, '  "Setup1 : Sweep", "Modal Solution Data", Array()\n');
% %end PCR
% %d_phase
% fprintf(fid, 'oModule.CreateOutputVariable "d_phase",  _\n');
% fprintf(fid, ' "ang_deg(S(FloquetPort1:2,FloquetPort1:2))-ang_deg(S(FloquetPort1:1,FloquetPort" & _\n');
% fprintf(fid, '"1:1))", "Setup1 : Sweep", "Modal Solution Data", Array()\n');
% %end
% %rll
% fprintf(fid, 'oModule.CreateOutputVariable "rll",  _\n');
% fprintf(fid, '  "1/2*(S(FloquetPort1:1,FloquetPort1:1)-S(FloquetPort1:2,FloquetPort1:2))-1j/2*(" & _\n');
% fprintf(fid, '  "S(FloquetPort1:1,FloquetPort1:2)+S(FloquetPort1:2,FloquetPort1:1))",  _\n');
% fprintf(fid, '  "Setup1 : Sweep", "Modal Solution Data", Array()\n');
% %end rll


%     hfssSaveProject(fid, PrjFile, true);
%     fprintf(fid, 'oDesign.AnalyzeAll\n');
%     导出数据
%      FileName = [tmpDataFile, 'VSWR', num2str(i)];
%      hfssExportToFile(fid, 'VSWR', FileName, Type);
% end
    fclose(fid);
    % 利用HFSS执行脚本程序
%     hfssExecuteScript(hfssExePath, tmpScriptFile);

% remove all the added paths.
 hfssRemovePaths(relPath);
 rmpath(relPath);
% add paths to the required m-files.

