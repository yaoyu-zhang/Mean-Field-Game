function [m0] = initialm0(delta_S,sigma)
%INITIALM0 此处显示有关此函数的摘要
%   此处显示详细说明
N = 1/delta_S;

m0 = 0:1:N;
m0 = m0/N;
m0 = 1/2/pi/(sigma^2)* exp(-(m0-0.5).^2/2/(sigma^2));
m0 = m0/sum(m0);

end

