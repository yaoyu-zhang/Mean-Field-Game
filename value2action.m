function [r] = value2action(v,p,rmin,rmax,delta_T,delta_S,gamma,Er)
rmin = rmin*delta_T/delta_S;
rmax = rmax*delta_T/delta_S;
[T,N] = size(v);
T = T-1;
N = N-1;
r = zeros(T,N+1);
for i = T:-1:1
    for j=1:1:N+1
        left = max(rmin,(1-j));
        right = min(rmax,(N+1-j));
        range = left:1:right;
        t = range*delta_S*Er/1000;
        cost1 = p(i)*( (t>0).*t/gamma+(t<=0).*t*gamma );
        cost2 = v(i+1,range+j);
        cost = cost1+cost2;
        [v(i,j),index] = min(cost);
        r(i,j) = range(index);
    end
end
r = r*delta_S/delta_T; 

end

