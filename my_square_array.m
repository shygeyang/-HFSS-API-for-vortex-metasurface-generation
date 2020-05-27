%目标：用于建立涡旋特性的metasurface.
%2020.5.26日 by Ling-Jun Yang

clear, clc

units = 'mm';
Units = 'mm';
Type = 'csv';

%%%%%%%需要的变量
fc=9;
M=8;
N=M;
arraynum=M*N;
%design Variables
    f=9;%GHz
    p=10;%单元周期
    h=3;
    b1=3;
    b2=8;
    s1=0.3;
    s2=3.5;
    da=12;
    ob_h=100;%观测高度
    
%%此段为计算不同模式阵列相位分布(squar_loop)。a矩阵为输出相位矩阵
clc , close all
a=ones(N,M);%a矩阵为输出相位矩阵

l=1.5;%想要的OAM模
Zf=M+M/2;
Zk=M/2;
Xk=Zk;
Yk=0;
k=1/2;
for m=1:M
    for n=1:N
      a(n,m)=l.*(angle(m-0.5-M/2+1j*(N/2-n+0.5))/pi*180);
      if a(n,m)<0
          a(n,m)=a(n,m)+360;
      end

    end
end
      a=k*a;
      a=rem(a,360);
      a=roundn(a,-2);
%%%%%%%%%脚本开始准备
false = 0;
true = 1;
fileName = ['squareloopN',num2str(N),'_',num2str(fc),'G_',num2str(arraynum),'N_l=',num2str(l)];
tmpScriptFile = ['D:\B站视频\超表面的建立\metasurface_HFSS_API\',fileName,'.vbs'];
%开始写VBS脚本
    fid = fopen(tmpScriptFile, 'wt');   % 'wt'表示以文本模式打开文件，可写，覆盖原有内容
   % 创建一个新的工程并插入一个新的设计
    hfssNewProject(fid);
    hfssInsertDesign(fid, fileName);
%     hfssPre(fid);
   
    hfssaddVar(fid, 'f0', f, []);
    hfssaddVar(fid, 'p', p, units);
    hfssaddVar(fid, 'h', h, units);
    hfssaddVar(fid, 'b1', b1, units);
    hfssaddVar(fid, 'b2', b2, units);  
    hfssaddVar(fid, 's1', s1, units);
    hfssaddVar(fid, 's2', s2, units);  
    hfssaddVar(fid, 'da', da, units); 
    hfssaddVar(fid, 'M', M, []);  
    hfssaddVar(fid, 'N', N, []); 
    hfssaddVar(fid, 'ob_h', ob_h, units); 
    %%%%%%%%Array of metasrface(square loop)
    squareloopName = cell(arraynum, 1);
    squareloopName1 = cell(arraynum, 1);
   cell_squareloopName=cell(1,1);
    for n=1:N
        for m=1:M
            i=N*(m-1)+n;
            %%%%add local CS
            CSName = ['CS', num2str(i)];
            Origin = {['(',num2str(m),' -(M+1)/2)*p'],['((N+1)/2-',num2str(n),')*p'],'0'};
            XAxisVec = [1 0 0 ];
            YAxisVec = {0, 1, 0};
            hfssCreateRelativeCS(fid, CSName, Origin,XAxisVec, YAxisVec, units);
            hfssSetWCS(fid, CSName);
   %%%%%%%Array of metasrface(square loop)
            squareloopName{i} = ['squareloop_',num2str(i)];
            Start_squareloop = {'-b1/2','-b2/2','h'};
            hfssRectangle(fid, squareloopName{i}, 'Z', Start_squareloop, 'b1', 'b2', units);

            squareloopName1{i} = ['squareloop1_',num2str(i)];
            Start_squareloop1 = {'-b1/2+s1','-b2/2+s2','h'};
            hfssRectangle(fid, squareloopName1{i}, 'Z', Start_squareloop1, 'b1-2*s1', 'b2-2*s2', units);
  
            hfssSubtract(fid, squareloopName{i}, squareloopName1{i}); 
            cell_squareloopName{1,1}= squareloopName{i};
            hfssRotate(fid,cell_squareloopName , 'Z', a(n,m));
            hfssSetWCS(fid, 'Global');
        end
    end
fclose(fid);
