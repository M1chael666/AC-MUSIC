%%%%%% Average Consensus函数（有限次数AC）%%%%%%%%%%%全联通
%%% 一共有P个节点，假设任意两个节点间都可以数据交换
%%% 数据存储在矩阵U中，节点之间的交换权值矩阵为W
%  输入：u0初始数据，行向量
%       节点个数P
%  输出：out AC值(经过Iac迭代后的数据求均值）
function out = AC(u0,P)

%%%%%% 先计算节点图的拉普拉斯矩阵
%%%对于一个连通图，任意两个节点之间都有边。生成P个节点连通图的边列表
s = [];
t = [];
for i = 1:P-1
    s = [s,i*ones(1,P-i)];
    t = [t,i+1:P];
end

G = graph(s,t);
L = laplacian(G);
[~,LAM] = eigs(L);
lambda = diag(LAM);

Rnum = min(size(LAM,1),P);
lambda = lambda(1:Rnum);

Iac = Rnum-1;
U = zeros(P,Iac);           %数据

%% 计算交换权值矩阵W        
lambda(lambda==0) = 1;
W1 = (-1)^(Rnum-1)/prod(lambda(1:Rnum-1))*(L-lambda(1)*eye(P));  %第一次迭代的权值

%%%第一次数据交换
U(:,1) = W1*u0.';

%%%Iac次数据交换
for i = 2:Iac
    U(:,i) = (L-lambda(i)*eye(P))*U(:,i-1);
end

out = sum(U(:,Iac));
end
