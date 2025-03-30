%%输入：
%   Pos_receive   初始位置
%   N_sample      采样数
%%输出：           移动后的坐标
%   P_receive_x   
%   P_receive_y   
%   P_receive_z
function [P_receive_x,P_receive_y,P_receive_z] = P_move(Pos_receive,N_sample)
%% 计算移动中的阵元位置
dread = pi/180;  %角度转弧度制
V_alpha =[196;145;208;133;138;-309;142;214;297;122;-148];   %速度的方向角
V = 7000;    %速度7000m/s
    x_dis = cos(V_alpha.*dread)*V/N_sample;
    y_dis = sin(V_alpha.*dread)*V/N_sample;
    z_dis = 0;
for ii = 1:N_sample
    P_receive_x(:,ii) = Pos_receive(:,1) + (ii-1)*x_dis; 
    P_receive_y(:,ii) = Pos_receive(:,2) + (ii-1)*y_dis; 
    P_receive_z(:,ii) = Pos_receive(:,3) + (ii-1)*z_dis; 
end
end


