%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%依次输入目标信号源的坐标（必须是单个信号源），接收阵列的坐标（阵元个数任意），信号波长
%输出该阵列对信号源的导向矢量(近场）
%注意theta是俯仰角度，phi是方位角度，坐标系是以参考接收阵元为原点(阵元1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = asteer_near(Pos_signal,Pos_receive,lambda)
%% 计算目标信号源的距离信息

R = (Pos_receive-(ones(size(Pos_receive,1),1)*Pos_signal)).^2;
R = sqrt(sum(R,2));
R_center = sqrt(sum(Pos_signal.^2,2));          %信号源到参考点（0，0，0）的距离
Rdiff = R-R_center;

out = exp(-1i*2*pi/lambda*Rdiff)/sqrt(size(Pos_receive,1));
end