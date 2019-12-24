clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        System identification (DS-CG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ensemble = 5;               % number of independent runs
N        = 10;                   % number of coefficients
K        = 2000;                      % number of iterations
gamma    = 10^(-4);                  % regularization factor
L        = 3;

lambda = 0.98;
eta = 0.48;
delta2 = 10^(-4);
sigma_ne = sqrt(1-0.88^2);
MSE = zeros(1088,11);
misalig_total = zeros(1991,11);
w_alll = zeros(8,5001);
sigma_k = sqrt(2);

uu = Mackey_glass;
uu = reshape(uu,1100,100);

for ee=1:ensemble

x = sqrt(0.001)*randn(1,K);
n = sqrt(0.001)*randn(1,K);

d(1) = -0.76*x(1) + 0.5*x(1)^2;
d(2) = -0.76*x(2) - 1*x(1) + 0.5*x(2)^2 - 1.6*x(1)^2;
for k = 3:K-9

    d(k) = -0.76*x(k)-x(k-1)+x(k-2)+0.5*x(k)^2+2*x(k)*x(k-2)-1.6*x(k-1)^2+1.2*x(k-2)^2+0.8*x(k-1)*x(k-2)+n(k);

end
num = 1;
den = [1 -0.95];
x   = filter(num,den,x);

n  = randn(size(uu,1),1);
u1 = uu(:,ee)+sqrt(0.01)*n;
x = u1(1:1100-L);
d = u1(N+L:end);

x_aux          =   x;
x_aux          =   buffer(x_aux,N,N-1,'nodelay');
x_tdl          =   flipud(x_aux);


cond1=0;
cond2=0;

P_up_all = 0.001:0.1*(1-0.001):1;

g=zeros(N,1);
c=zeros(N,1);

for p = 1:length(P_up_all)
    
P_up =  P_up_all(p);                
r    =  zeros(1088,2000); 
beta = ones(1,1991);
G = km_kernel(x_tdl,x_tdl,'gauss',sigma_k);
r(:,1)    = d;
r_old     = r(:,1);
v    = G*r(:,1);
u       = zeros(1991,1);
beta(1) = 1;
eta  = r(:,1)'*G*r(:,1); 

threshold_tau = floor(0.02*size(x_aux,2));
teste1 = [];
e_antigo = 0;
 b = 0.99;
 
for k=1:size(x_aux,2)

tau_max = 40000;
  if  k < threshold_tau
      sqrt_tau_max1 = inf;
      sqrt_tau_max2 = -inf;
  end
  if k >= threshold_tau
      sqrt_tau_max1 = mean(teste1) + 3*sqrt(var(teste1));
      sqrt_tau_max2 = mean(teste1) - 3*sqrt(var(teste1));
  end  
sqrt_tau_max1 = sqrt(40000000000000000000000);
sqrt_tau_max2 = -sqrt(40000000000000000000000);

   alpha2 = P_up*N*(1-lambda)/(2-P_up*(1-lambda));

 sqrt_tau = sqrt(1+alpha2)*qfuncinv((P_up)/2);
   sqrt_tau = sqrt(1+alpha2)*qfuncinv((P_up+2*qfunc(sqrt_tau_max1/sqrt(1+alpha2)))/2);

    e(k) = d(k) - km_kernel(x_tdl(:,k),x_tdl,'gauss',sigma_k)*r(:,1:k)*u(1:k);
    e_atual = e(k)^2;
    teste(k) = abs(e(k))/sqrt((1-b)*e_atual+b*e_antigo);
    teste2(k) = (e(k))/sqrt(0.001);
    
    if (teste(k) <= sqrt_tau)
        delta(k)=0;
        cond1=cond1+1;
        
    elseif (teste2(k) > sqrt_tau_max1) || teste2(k) < sqrt_tau_max2
        delta(k)=0;
        cond2 = cond2+1;
        
    else
        delta(k)=1;
    end
   
    if (delta(k) == 0)
        alpha(k) = 0;
        teste1 = [teste1 e(k)/sqrt(0.001)];    
        if (teste2(k) > sqrt_tau_max1) || teste2(k) < sqrt_tau_max2  
            e(k) = 0.0;
            d(k) = 0;
        end
    else
    
    teste1 = [teste1 e(k)/sqrt(0.001)];    
    
    alpha(k) = eta/dot(v,v);
    pi       = 1;
    
    for j=k:-1:1
        u(j) = u(j) + alpha(k)*pi;
        pi   = pi*beta(j);
    end
      
    r(:,k+1)    = r_old - alpha(k)*v;
    eta_old = eta;
    eta  = dot(r(:,k+1),G*r(:,k+1));
    beta(k+1) = eta/eta_old;
    v    = beta(k+1)*v + G*r(:,k+1);
    r_old = r(:,k+1);
    
    end
     e_antigo = (1-b)*e_atual+b*e_antigo;

end
MSE(:,p)= MSE(:,p)+(abs(e.')).^2;
P_up_est(p) = size(delta(delta ==1),2)/(size(delta,2));
e_all(:,p) = e;
tau(p) = sqrt_tau;
d_all(:,p) = d;
dd_all(:,p) = di;
y_all(:,p) = y;
if p == 6
    w_all = w;
end

end
w_alll = w_alll + w_all;
P_uppp(:,ee)=P_up_est;
ee
end

P_uppp = mean(P_uppp');
MSE = MSE/ensemble;
misalig_total = misalig_total/ensemble;

figure,
plot(10*log10(MSE(:,1)));
set(gca,'fontsize',18)
xl = xlabel('Number of iterations, $k$'); 
yl = ylabel('Misalignment (dB)'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
saveas(gcf,'id_data_set2','fig');
