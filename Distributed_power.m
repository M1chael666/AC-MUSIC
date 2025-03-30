%%%%%%%分布式的功率计算%%%%%%旧
%%%在每个节点中分布式估计协方差矩阵
%  输入：X  阵列输出
%        P  节点个数
%        K  采样快拍数
%        M  阵元个数
%        Ipm 迭代次数
%  输出：out Eg_1 特征向量迭代终值
%
function out = Distributed_power(X,P,K,M,Ipm)
Eg = zeros(M,M);          %第n次的值存为Eg
Eg_1 = zeros(M,M);        %第n+1次的值存为Eg_1
fi = zeros(1,M);

for ii = 1:Ipm
    if ii == 1
        Eg = randn(M,M);
    else
        Eg = Eg_1;        %数据更新
    end                   %%%%%%%%%第一步迭代随机取初值，后续初值取上一步
    
    
    Eg_1 = zeros(M,M);    %Eg_1归零
    for m = 1:M 
       fi(m) = AC(Eg_1(:,m)'*diag(f11(X,P,K,M,Eg(:,m))),P);   %存疑
       Eg_1(:,m) = f11(X,P,K,M,Eg(:,m))-P*sum(Eg_1*diag(fi),2);   %公式16
       
       if ii == Ipm
           Eg_1(:,m) = Eg_1(:,m)/abs(P*AC(Eg_1(:,m)'*diag(Eg_1(:,m)),P))^0.5;     %归一化
       else
       end
    end
    
    out = Eg_1;
end