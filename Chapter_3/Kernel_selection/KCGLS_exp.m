function [ e] = KCGLS_exp( x,d,S )
N = S.filterOrderNo;
K = S.iterationsnumbers;
lambda = 0.99;

x_aux          =   x;
x_aux          =   buffer(x_aux,N,N-1,'nodelay');
x_tdl          =   flipud(x_aux);
p = 0.0001;

w(:,1)    = zeros(N,1);
g    = zeros(N,1);%X*r;
c    = g;
% R = zeros(N); 
eta = 0.49; 

for k=1:1000
    e(k)      = (d(k) - w(:,k)'*x_tdl(:,k));
%      R = lambda*R+x_tdl(:,k)*x_tdl(:,k)';
    if k==1
        X  = [x_tdl(:,k)];
    else 
        X  = [sqrt(lambda)*X x_tdl(:,k)];
    end
    v = X'*c;
    alpha(k) = eta*(c'*g)/(v'*v+p);  
    g1 = g;
    
    if k==1
        r = -alpha(k)*v;
    else
        r = lambda*[r;0]-alpha(k)*v;
    end
    g = X*r+x_tdl(:,k)*(e(k));

%     g = lambda*g-alpha(k)*(X*X')*c+x_tdl(:,k)*(e(k));%
    w(:,k+1) = w(:,k) + alpha(k)*c;
    beta(k) = ((g-g1)'*g)/(g1'*g1+p);
    c  = g + beta(k)*c; 

end