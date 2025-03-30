%%%%%%分布式的music算法
%%%分布式的music算法的区别在于不再计算整体的协方差矩阵Rxx
%%%省去了将各阵元的输出传输到某一中心的过程
%%%这里将单个阵元看作一个节点
%%%参考文献Distributed_Algorithms_for_Array_Signal_Processing
%%%公式（9）（11）不需要整体的Rxx，只需要每一时刻的节点输出就可以得到准确的特征向量
%   输入：
%
%
%
%
%
%   输出：
%
function [X,Y,Z] = Distributed_music(Pos_signal,Pos_receive,X_search,Y_search,lambda,snr,fig_mark)
K = size(Pos_signal,1);                        %信号源个数                                           %信源数
M = size(Pos_receive,1);                       %阵元个数                                           %阵元个数
P = M;                                         %节点个数
Z_search = Pos_signal(1,3);                    %信号源Z坐标


%% search arrival direction
N_sample = 10;                                 %采样点数
signal = sig_generation(N_sample);             %LFM信号                                                             
                                                
A_MxK = zeros(M,K);                        %近场导向矢量  
for k = 1:K
    A_MxK(:,k) = asteer_near(Pos_signal(k,:),Pos_receive,lambda);
end                                      


%%%%%%%%%%   信号接收
c = 3e8;
f0 = c/lambda;
yin = A_MxK*signal;                     
% P_sig = trace(yin*yin');
yin = awgn(yin,snr,'measured');         %加入噪声

Ipm = 20;                  %迭代次数
Un = Distributed_power_1(yin,P,N_sample,M,Ipm);  %使用分布式功率计算噪声特征矢量

P_MUSIC = zeros(length(X_search),length(Y_search));
for ii = 1:length(X_search)
    for jj = 1:length(Y_search)
        asteer = asteer_near([X_search(ii),Y_search(jj),Z_search],Pos_receive,lambda);
        Pow = 0;
        for k = 1:K
            Pow = Pow+abs(P*AC(Un(:,k)'*diag(asteer),P))^2;
        end
        P_MUSIC(ii,jj) = 1/M-Pow;          %公式33
    end
end


if(fig_mark)
figure
mesh(Y_search,X_search,abs(P_MUSIC));
title('MUSIC谱','FontSize',13);
xlabel('Y方向/m','FontSize',13);
ylabel('X方向/m','FontSize',13);
zlabel('空间功率/dB')
else
end


%%%%%%%%%%% 搜索谱峰
[max1,locs1] = max(abs(P_MUSIC));
[~,locs2] = max(max1);
Y = Y_search(locs2);
X = X_search(locs1(locs2));
Z = Z_search;
end