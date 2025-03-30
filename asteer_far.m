%%计算目标信号源的距离信息
% Rm_s_est(i,j) = S_R- (A_X(i,j)*S_X+A_Y(i,j)*S_Y+A_Z(i,j)*S_Z)/S_R;
%A(i,j) = exp(-1i*2*pi/w_l*(Rm_s_est(i,j)-S_R));

function out  = asteer_far(Pos_signal,Pos_receive,lambda)
u = sig_u(Pos_signal);      %单位向量 列向量
out  = exp(1i*2*pi/lambda*Pos_receive*u);
end
