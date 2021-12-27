function [r] = value2action(v,p,rmin,rmax,delta_T,delta_S,lambda,S_ev,gamma)
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
        t = range*delta_S*S_ev/1000;
        cost0 = -1/2*gamma*(t.^2)+S_ev*t;
        cost1 = -p(i)*( t + t.^2*lambda);
        cost2 = v(i+1,range+j);
        cost = cost0+cost1+cost2;
        [v(i,j),index] = max(cost);
        r(i,j) = range(index);
    end
end
r = r*delta_S/delta_T; 

end

