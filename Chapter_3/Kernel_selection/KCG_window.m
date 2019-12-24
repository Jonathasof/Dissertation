function [e] = KCG_window( x,d,S )
N = S.filterOrderNo;
M = 5;
x1 = zeros(N,M-1);
delta2 = 10^(-4);

x_aux          =   x;
x_aux          =   buffer(x_aux,N,N-1);
x_tdl          =   [x1,flipud(x_aux)];

G       = km_kernel(x_tdl(:,1:M),x_tdl(:,1:M),'gauss',S.sigma_k);
G(1:end-1,:)       = 0;
G(end,1:end-1)     = 0;
r(:,1)  = [zeros(1,M-1), d(1)];
u       = zeros(S.I_max,1);
beta(1) = 1;
eta(1)  = r(:,1)'*G*r(:,1); 
G       = km_kernel(x_tdl(:,1+1:1+M),x_tdl(:,1+1:1+M),'gauss',S.sigma_k);
G(1:end-2,:)       = 0;
G(:,1:end-2)       = 0;
v(:,1)  = G*r(:,1);

for i=1:S.I_max-2
    alpha(i) = dot(r(:,i),G*r(:,i))/dot(v(:,i),v(:,i));%0.5%
    pi       = 1;
    
    for j=i:-1:1
        u(j) = u(j) + alpha(i)*pi;
        pi   = pi*beta(j);
    end
    
    r(:,i+1)  = r(:,i) - alpha(i)*v(:,i);
    eta(i+1)  = dot(r(:,i+1),G*r(:,i+1));
    beta(i+1) = (eta(i+1))/(eta(i));
    %-dot(r(:,i),km_kernel(x_tdl(:,i+1:i+M),x_tdl(:,i+2:i+M+1),'gauss',S.sigma_k)*r(:,i+1))
    if i<=M-3
        for j = i:-1:1
            w    = 0;
            w    = w + km_kernel(x_tdl(:,i+M),x_tdl(:,j+1:j+M),'gauss',S.sigma_k)*r(:,j)*u(j);
        end
%         w    = (km_kernel(x_tdl(:,i+M),x_tdl(:,i+1:i+M),'gauss',S.sigma_k))*r(:,1:i)*u(1:i);
        e(i) = d(i+1) - w;
        G    = km_kernel(x_tdl(:,i+2:i+M+1),x_tdl(:,i+2:i+M+1),'gauss',S.sigma_k);
        G(1:end-(i+2),:)       = 0;
        G(:,1:end-(i+2))       = 0;
    else
        for j = i:-1:1
            w    = 0;
            w    = w + km_kernel(x_tdl(:,i+M),x_tdl(:,j+1:j+M),'gauss',S.sigma_k)*r(:,j)*u(j);
        end
%         w    = km_kernel(x_tdl(:,i+M),x_tdl(:,i+1:i+M),'gauss',S.sigma_k)*r(:,1:i)*u(1:i);
        e(i) = d(i+1) - w;
        G = km_kernel(x_tdl(:,i+2:i+M+1),x_tdl(:,i+2:i+M+1),'gauss',S.sigma_k);
    end
    v(:,i+1) = beta(i+1)*v(:,i) + G*r(:,i+1);

end
