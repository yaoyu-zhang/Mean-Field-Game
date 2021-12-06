load ERCOT; % all data
close all;
Di = demand; % Load demand data other than EV
delta_S = 0.0025; 
Time_length = length(Di); 
delta_T = 1;
T = Time_length/delta_T;
N = 1/delta_S; 
new_Di = zeros(T,1); 
for i =1:1:length(Di)
    for j=1:1:1/delta_T
        new_Di(i/delta_T-1/delta_T+j) = Di(i)*delta_T;
    end
end
Di = new_Di; 
S_all = 25; 
P_ev = 2.5 ; S_ev=25;
rmin = -P_ev/S_ev;
rmax = P_ev/S_ev;
gamma = 0.95; 

sigma = 1.2;
m0 = initialm0(delta_S,sigma); 
m = ones(T,1)*m0;
new_m = m;
c = 10000000;
phi = c*( (0:1:N)/N-1).^2; %ÐèÒª¸ø¶¨ 1ÐÐ£¬N+1ÁÐ
V = zeros(T+1,N+1);
V(T+1,:) = phi;
r = zeros(T,N+1);

cul = [];
ccul = [];
epsilon_p = 0.00001;
epsilon_D = 0.00001;
epsilon = 0.00001;
new_Ds = zeros(T,1);
Ds = new_Ds + 2*epsilon_D;
circle1 = 1;
kk = 1;
while exp(-0.01*kk)*abs(Ds-new_Ds)> epsilon_D
    tic;
    Ds = (1-exp(-0.01*kk) )*Ds + exp(-0.01*kk)*new_Ds;
    kk = kk+1;
    m = new_m;
    new_p = f_price_ERCOT(Ds+Di,C,MF,delta_T); 
    p = new_p + 2*epsilon_p;
    cul = [];
    circle2 = 1;
    while sum(abs(p-new_p)) > epsilon_p
        
        p = new_p;

        r = value2action(V,p,rmin,rmax,delta_T,delta_S,gamma,S_ev);
        new_p = f_price_ERCOT(  sum(m.*( (r>0).*r/gamma+ (r<=0).*r*gamma  ),2)*delta_T*S_all    +Di,C,MF,delta_T);
        cul = [cul,sum(abs(p-new_p))];
    end
    new_m = action2m(r,m0,delta_T,delta_S);
    new_Ds = sum(new_m.*((r>0).*r/gamma+ (r<=0).*r*gamma),2)*delta_T*S_all;
    ccul = [ccul,sum(abs(Ds-new_Ds))];
    toc;
end
m =new_m;
Ds = new_Ds;
new_p = f_price_ERCOT(  sum(m.*( (r>0).*r/gamma+ (r<=0).*r*gamma  ),2)*delta_T*S_all    +Di,C,MF,delta_T);
%-------------------------------
profit = zeros(1,N+1);
action = zeros(T,N+1);
state = zeros(T+1,N+1);
for j=1:1:N+1
    state(1,j) = j;
    for i = 2:1:(T+1)
        action(i-1,j) = r(i-1,state(i-1,j));
        state(i,j) = state(i-1,j)+ round(action(i-1,j) *delta_T/delta_S);
    end
    profit(j) = sum(new_p.*action(:,j));
end
action = ((action>0).*action/gamma+ (action<=0).*action*gamma)*delta_T*S_ev/1000;
for j= 1:1:N+1
    profit(j) = sum(new_p.*action(:,j));
end

 save("EV_mfg");

