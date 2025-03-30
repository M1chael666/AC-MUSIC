%%输入    x,y,z坐标 (默认为一个信号源)
%%输出    单位向量u 列向量
function out = sig_u(Pos_signal)
   R = sqrt(sum(Pos_signal.^2,2));
   out = Pos_signal'/R;
end