 %%%%%%公式11的计算函数 参考文献Distributed_Algorithms_for_Array_Signal_Processing
      %   输入：X 阵列输出
      %         P 节点个数
      %         K 快拍数
      %         M 阵元个数
      %         e_n 迭代对象，列向量
      %   输出：result 下一步的迭代结果
function out = f11(X,P,K,M,e_n)
     b = zeros(1,K);
     result = zeros(M,K);
     for k = 1:K
         b(k) = AC_1(X(:,k)'*diag(e_n),P);
              %%%这里需要整体计算AC          
         result(:,k) = X(:,k)*b(k);   %这一步是在每个节点中分别计算的
     end

     out = sum(result,2)*P/K;   %这一步是在每个节点中分别计算的

end
