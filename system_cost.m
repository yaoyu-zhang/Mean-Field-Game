function [cost] = system_cost(C,MF,D,delta_T)
%SYSTEM_COST 此处显示有关此函数的摘要
%   此处显示详细说明
cost = zeros(size(D));
D = D/delta_T;
for j = 1:1:length(D)
    temp = zeros(size(C));
    for i = 1:1:length(C)
        if i==1
            temp(i) = C(i);
        else
            temp(i) = C(i)-C(i-1);
        end
        
        if D(j)>C(i) && D(j)<=C(i+1)  
            temp(i+1) = D(j) - C(i);
            break;
        end
    end
    cost(j) = sum(MF.*temp*1000);
end
cost = sum(cost);
end

