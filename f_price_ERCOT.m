function [p] = f_price_ERCOT(D,C,MF,delta_T)
% 同意采取GWh为单位
p = zeros(size(D));
D = D/delta_T;
% D = D./delta_T;
for j = 1:1:length(D)
    for i = 1:1:length(C)
        if D(j)>C(i) && D(j)<=C(i+1)  
            break;
        end
    end
    p(j) = MF(i);
end
% p = MF(i);
%F_PRICE 此处显示有关此函数的摘要
%   根据需求返回价格
% p = 2*D.^2+20;
% p = D;
% p = p;
end

