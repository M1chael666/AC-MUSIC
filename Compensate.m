%%输入：
%导向矢量：


%%输出：
%补偿后导向矢量：
function out = Compensate(Pos_signal,A_MxK,N_sample,lambda)
dread = pi/180;  %角度转弧度制
V_alpha =[196;145;208;133;138;-309;142;214;297;122;-148];   %速度的方向角
V = 7000;   %阵元速度
t = 0.01;   %相邻采样时间间隔
x_dis = cos(V_alpha.*dread)*V*t;        %x方向位移
y_dis = sin(V_alpha.*dread)*V*t;        %y方向位移
z_dis = zeros(length(x_dis),1);         %z方向位移
R_dis = [x_dis,y_dis,z_dis];            %相邻采样时刻阵元位移

u = sig_u(Pos_signal);                  %单位向量 列向量
com = exp(1i*2*pi/lambda*R_dis*u);      %相位补偿因子   列向量

for ii = 1:N_sample
    A_MxK_com(:,ii) = A_MxK(:,ii).*(com.^(N_sample-ii));    %迭代补偿
end

out = A_MxK_com;
end