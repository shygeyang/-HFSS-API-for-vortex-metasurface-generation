%CLACULATE_FIELD 此处显示有关此函数的摘要
%   此处显示详细说明         $$$构造一定donut形状的超材料  幅度相位调控matlab计算场
clc
clear all
close all
R0=0.06;%建模的中心半径  单位m
R1=0.06;%建模的波动半径   单位m
d=0.01;  %周期大小   单位m
rn=(R0+R1)/d;%建模半径rn个面积内建模R0+R1个周期
l=[11];%想要的OAM模,可以为数组
mag_l=[1];%想要的OAM模对应的模式比重,可以为数组
% distance,  距离   单位m
f=19*10^9; %  频率  Hz
% [sample_m],  采样阵列大小个数
% [sp]  采样周期大小  单位m
fto=[15];%GHz

for fi=1:length(fto)

  f=fto(fi)*10^9; %  频率  Hz  
%准备工作
theta=0:2*pi/1000:2*pi;
a1_total=zeros(1,length(theta));%需要的复数OAM场
for o=1:length(l)
    a1_total= a1_total+mag_l(o)*exp(-1j*l(o)*theta);
end
max_a= max(abs(a1_total));%需要的复数OAM场的最大值用于归一化
% cin=12-4/max_a.*abs(a1_total)-0.5;
% cin=12+4/max_a.*abs(a1_total)+0.5;
% phase=angle(al_total);

N=ceil(rn*2/sqrt(3));%最大建模N圈六边形
% N=10;
% arraynum=6*(1+N)*N/2;%最多单元个数

k=2*pi*f/(3*10^8);

%求解矩阵sampling_m
theta_S=0:1:90;
phi_S=1:1:360;

sampling_matrix=zeros(length(theta_S),length(phi_S));
% phi0=phi0/180*pi;

for o=1:1:length(theta_S)
    for p=1:1:length(phi_S)
   
        for n=1:N %圈数
            for j=1:6 %边
               for i=1:n %个
                     x=n*d*cos((j-1)*pi/3)-(i-1)*cos((j-2)*pi/3)*d;%%描述n,j,i,时刻x，y轴的坐标
                     y=n*d*sin((j-1)*pi/3)+(i-1)*sin((j+1)*pi/3)*d;
                     al=angle(x+1j*y)/pi*180;%%%描述n,j,i,时刻方位角al
                    if al<0
                              al=al+360;
                    end
                 a_total=0;%累加需要的幅度和相位A总
                    for m=1:length(l)
                        a_total= a_total+mag_l(m)*exp(-1j*l(m)*al*pi/180);%该方位角al的复数值a_total。
                    end
                a_total=a_total/max_a;%调整幅度归一化
                a_=angle(a_total);%该单元需要的相位rad

                    %判断是否需要建模
                     if sqrt(x^2+y^2)<abs(a_total)*R1+R0 && sqrt(x^2+y^2)>R0-abs(a_total)*R1     %abs(a_total)

%                          r=sqrt((x+(o-0.5-sample_m/2)*sp)^2+(y+(sample_m/2-p+0.5)*sp)^2+distance^2);
%                             r=sin(theta_S(o)/180*pi)*cos(phi_S(p)/180*pi)*x+sin(theta_S(o)/180*pi)*sin(phi_S(p)/180*pi)*y;
%                          sampling_matrix(o,p)=(1+cos(theta_S(o)/180*pi))*exp(-1j*(k.*r)-1j*a_)+sampling_matrix(o,p);
                           r=sin(theta_S(o)/180*pi)*cos(phi_S(p)/180*pi)*x+sin(theta_S(o)/180*pi)*sin(phi_S(p)/180*pi)*y;
                         sampling_matrix(o,p)=(1+cos(theta_S(o)/180*pi))*exp(-1j*(k.*r)-1j*a_)+sampling_matrix(o,p);
                     end
               end
            end
        end      
     end
end
mag_s=abs(sampling_matrix);
angle_s=angle(sampling_matrix);
[A,B]=max(mag_s);
thetaM=mean(B)-1;

% [X,Y] = meshgrid(theta_S,phi_S);
[X,Y] = meshgrid(phi_S,fliplr(theta_S));
figure(1)
mesh(X,Y,mag_s)
colormap(hot)
 colorbar
 view(0,90)
figure(2)
mesh(X,Y,angle_s)
load mycolor
colormap(rad2blue)
view(0,90)
 colorbar



end

