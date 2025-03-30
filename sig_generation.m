%%%%%%%%%%目标辐射源发射的信号
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入快拍数，信噪比
%输出采样之后的信号波形
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sig_generation(N_sample)
%%   产生LFM标准发射信号 线性调频信号

% 载频：f0=0;
% 脉宽：pulsewidth=200;微秒
% 频率间隔：fdelta=0.5;M
% 带宽：bandwidth=0.5;M
% 采样频率：fs=20;M
% 输出：单个通道的标准信号S(t)
bandwidth = 0.5;                     
pulsewidth = 20000;                  
fs = 20;
FB = bandwidth/2 ;                   
j = sqrt(-1);                        
miu = bandwidth/pulsewidth;          %频率变化率
time = 0:1/fs:pulsewidth-1/fs;

signal = exp(j*2*pi*(0.5*miu*time.^2-FB*time));

out = signal(1:N_sample);
end
