function [ee,d_til] = KCG_exponential_2( x,d,S )
N = S.filterOrderNo;
K = S.iterationsnumbers;
lambda = sqrt(1);

x_aux          =   x;
x_aux          =   buffer(x_aux,N,N-1);
x_tdl          =   flipud(x_aux);

x_test = x_tdl(:,1001:end);
x_tdl  = x_tdl(:,1:1000);

d_test = d(1001:end);
d      = d(1:1000);

type = 'gauss';
X       = x_tdl(:,1);
q{1}  = km_kernel(X,X,type,S.sigma_k);
G_t = q{1};
G_M = G_t;
eta{1}  = d(1)/q{1};
d_til(1) = d(1);
e{1} = 0;
G1      = 1; 
   
for k=1:990
    clear r
    clear v
    g_k  = km_kernel(X,x_tdl(:,k+1),type,S.sigma_k);
    q_k  = km_kernel(x_tdl(:,k+1),x_tdl(:,k+1),type,S.sigma_k);
    d_til(k+1) =  km_kernel(X,x_tdl(:,k+1),type,S.sigma_k)'*(eta{k});
    ee(k) = (d(k+1)' - d_til(k+1)).^2;%%%%
    e{k+1}  = d(1:k+1)' - (G1.*km_kernel(X,x_tdl(:,1:k+1),type,S.sigma_k))'*eta{k}; %%%%
    X = [X x_tdl(:,k+1)];
    G_t     = [G_t g_k ; g_k' q_k];
%     r(:,1) = [ e{k}; d(k+1) - (G1.*g_k')*eta{k}];%%%%%%%%%%5
    G1    = [lambda*G1;1];
    G_M   = (G1*G1').*G_t;
    r(:,1) = e{k+1};
    gamma(1) = dot(G_M*r(:,1),r(:,1));
    v(:,1)   = G_M*r(:,1);
    beta(1) = 1;
    u       = zeros(2,1);
    for j = 1:2
        alpha(j) = gamma(j)/dot(v(:,j),v(:,j));
        pi       = 1;
    
        for l=j:-1:1
            u(l) = u(l) + alpha(j)*pi;
            pi   = pi*beta(l);
        end
        
        r(:,j+1) = r(:,j) - alpha(j)*v(:,j);
        gamma(j+1) = dot(r(:,j+1),G_M*r(:,j+1));
        gamma1     =  dot(r(:,j),G_M*r(:,j+1));
        beta(j+1)  = (gamma(j+1))/(gamma(j));
        v(:,j+1)   = (beta(j+1)*v(:,j)) + G_M*r(:,j+1);
        
    end
    
    eta{k+1} = [eta{k}; 0] +  (G1).*((r(:,1:end-1))*u);%((G1)).*
%     e{k+1}   = d(1:k+1)' - (G_M)*eta{k+1};
%     ker_D_X = km_kernel(X, x_test,'gauss',S.sigma_k);
%     d_hat = ((eta{k+1})'*ker_D_X);
%     err = d_test - d_hat;
%     ee(k) = mean(err.^2);
%     ee(k) = (d(k+1) - km_kernel(X,x_tdl(:,k+1),'gauss',S.sigma_k)'*((G1).*eta{k}))^2;
end

end
