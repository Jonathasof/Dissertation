function [ e] = KCG( x,d,S )
N = S.filterOrderNo;
K = S.iterationsnumbers;

x_aux          =   x;
x_aux          =   buffer(x_aux,N,N-1);
x_tdl          =   flipud(x_aux);
S.I_max = 1000;
G = km_kernel(x_tdl,x_tdl,'linear',S.sigma_k);
r(:,1)    = d;
v(:,1)    = G*r(:,1);
u       = zeros(S.I_max,1);
beta(1) = 1;
eta(1)  = r(:,1)'*G*r(:,1); 

for i=1:S.I_max
    alpha(i) = eta(i)/dot(v(:,i),v(:,i));
    pi       = 1;
    
    for j=i:-1:1
        u(j) = u(j) + alpha(i)*pi;
        pi   = pi*beta(j);
    end
    
    r(:,i+1)    = r(:,i) - alpha(i)*v(:,i);
    eta(i+1)  = dot(r(:,i+1),G*r(:,i+1));
    beta(i+1) = eta(i+1)/eta(i);
    v(:,i+1)    = beta(i+1)*v(:,i) + G*r(:,i+1);
    e(i)      = d(i) - km_kernel(x_tdl(:,i),x_tdl,'linear',S.sigma_k)*r(:,1:i)*u(1:i);
end
