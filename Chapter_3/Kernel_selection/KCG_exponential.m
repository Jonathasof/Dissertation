function [e] = KCG_exponential( x,d,S )
N = S.filterOrderNo;
lambda = sqrt(1);
delta2 = 10^(-4);

x_aux          =   x;
x_aux          =   buffer(x_aux,N,N-1);
x_tdl          =   flipud(x_aux);
[~,M]          =   size(x_tdl);

g       = [zeros(N,M-1), x_tdl(:,1)];
G       = zeros(M);
G(M,M)  = 1;
% G  = km_kernel(g,g,'gauss',S.sigma_k);
r(:,1)  = [zeros(1,M-1) d(1)];%%%
u       = zeros(S.I_max,1);
beta(1) = 1;
eta(1)  = r(:,1)'*G*r(:,1);  
gg{1}   = [g(:,2:end), x_tdl(:,2)];
G       = km_kernel(gg{1},gg{1},'gauss',S.sigma_k);
% G       = [zeros(M-2,M-2) zeros(M-2,2); zeros(2,M-2) km_kernel(x_tdl(:,1:2),x_tdl(:,1:2),'gauss',S.sigma_k)];
[g1,~]  = size(G);
G1      = [lambda*ones(1,g1-1),1]; 
G_t     = (G1'*G1).*G; 
v(:,1)  = G_t*r(:,1);
lamb(1,:)    = [lambda*ones(1,g1-1), 1];
w       = 0; 

for i=1:S.I_max-10
    e(i)     = d(i+1) - w;
    alpha(i) = eta(i)/(dot(v(:,i),v(:,i)));
    pi       = 1;
    
    for j=i:-1:1
        u(j) = u(j) + alpha(i)*pi;
        pi   = pi*beta(j);
    end
    
    r(:,i+1)  = r(:,i) - alpha(i)*v(:,i);
    eta(i+1)  = dot(r(:,i+1),G_t*r(:,i+1));
    beta(i+1) = eta(i+1)/(eta(i));
%     e(i)      = d(i+1) - [zeros(1,1088-(i+1)) ones(1,(i+1))].*[lamb.*km_kernel(x_tdl(:,i+1),g,'gauss',S.sigma_k)]*r(:,1:i)*u(1:i);%%%

    for j=i:-1:1
        w         = 0;
        w         = w + lamb(j).*km_kernel(x_tdl(:,1+i),gg{j},'gauss',S.sigma_k)*r(:,j)*u(j);
    end

%     w         = [lamb(i,:).*km_kernel(x_tdl(:,i+1),gg{i},'gauss',S.sigma_k)]*r(:,1:i)*u(1:i); %[zeros(1,1088-(i+1)) ones(1,(i+1))].*
    llamb     = lamb(i,:);
    lamb(i+1,:)      = [lambda*llamb(2:end), 1];
    gg{i+1}         = [gg{i}(:,2:end), x_tdl(:,i+2)];
    temp      = fliplr(km_kernel(x_tdl(:,i+2),gg{i+1},'gauss',S.sigma_k)); 
    G         = [G(2:end,2:end) temp(1:end-1)'; temp];
    G1        = [lambda*G1(2:end),1]; 
    G_t       = (G1'*G1).*G; 
    v(:,i+1)  = beta(i+1)*v(:,i) + G_t*r(:,i+1);
end

end
