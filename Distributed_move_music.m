%%分布式动平台MUSIC算法

function [X,Y,Z] = Distributed_move_music(Pos_signal,Pos_receive,X_search,Y_search,lambda,snr,fig_mark)
K = size(Pos_signal,1);                        %信号源个数                                           %信源数
M = size(Pos_receive,1);                       %阵元个数                                           %阵元个数
P = M;                                         %节点个数
Z_search = Pos_signal(1,3);                    %信号源Z坐标


%% search arrival direction                                     
N_sample = 100;                               %每秒钟采样点数                                     
signal = sig_generation(N_sample);             %LFM信号                                                             

[receive_x,receive_y,receive_z] = P_move(Pos_receive,N_sample);
A_MxK_move = zeros(M,N_sample);            %导向矢量  
for n = 1:N_sample                         
    A_MxK_move(:,n) = asteer_far(Pos_signal,[receive_x(:,n),receive_y(:,n),receive_z(:,n)],lambda);
end
   
                                     
%%%%%%%%%%   信号接收
c = 3e8;
f0 = c/lambda;
       

yin = A_MxK_move*diag(signal);                 %接收信号
yin = awgn(yin,snr,'measured');                %加噪声
Ipm = 10;                                       %幂迭代次数   


yin_com  = Compensate(Pos_signal,yin,N_sample,lambda);
Un = Distributed_power_1(yin_com,P,N_sample,M,Ipm);
         

P_MUSIC = zeros(length(X_search),length(Y_search));
for ii = 1:length(X_search)
    for jj = 1:length(Y_search)
        %yin_com  = Compensate([X_search(ii),Y_search(jj),Z_search],yin,N_sample,lambda);
        %Un = Distributed_power_1(yin_com,P,N_sample,M,Ipm);
        asteer = asteer_far([X_search(ii),Y_search(jj),Z_search],[receive_x(:,N_sample),receive_y(:,N_sample),receive_z(:,N_sample)],lambda);
        Pow = 0;
        for k = 1:K
            Pow = Pow+abs(P*AC_1(Un(:,k)'*diag(asteer),P))^2;
        end
        P_MUSIC(ii,jj) = 1/M-Pow;          %公式33
    end
end


if(fig_mark)
figure
mesh(Y_search,X_search,abs(P_MUSIC));
xlabel('y/(m)','FontSize',13);
ylabel('x/(m)','FontSize',13);
zlabel('谱功率/dB')
else
end


%%%%%%%%%%% 搜索谱峰
[max1,locs1] = max(abs(P_MUSIC));
[~,locs2] = max(max1);
Y = Y_search(locs2);
X = X_search(locs1(locs2));
Z = Z_search;
end