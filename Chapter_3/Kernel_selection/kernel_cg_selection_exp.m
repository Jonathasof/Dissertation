clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         System identification (DS-CG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('data.mat')
% data_x = data_wind;
% data_x = (data_x - mean(data_x))/var(reshape(data_x,65536*5,1))^(1/2);
    
ensemble = 5;               % number of independent runs
N        = 10;                   % number of coefficients
K        = 2000;                      % number of iterations
L        = 3;

lambda = sqrt(1);
eta1 = 1;
MSE = zeros(1000,11);
% misalig_total = zeros(1991,11);
w_alll = zeros(8,5001);
sigma_k = 1;

uu = Mackey_glass;
uu = uu + sqrt(0.001)*randn(size(uu));

for ee=1:ensemble
    
%% identificatxion 
x = sqrt(0.1)*rand(1,1100);
n = randn(1,1100);

den = 1;
num = [1 -0.5];% 0.3 0.7 -0.3];
x1   = filter(num,den,x); %Filter 

d(1) = -0.76*x1(1) + 0.5*x1(1)^2; %nonlinear mapping
d(2) = -0.76*x1(2) - 1*x1(1) + 0.5*x1(2)^2 - 1.6*x1(1)^2;
for k = 3:1100

    d(k) = -0.76*x1(k)-x1(k-1)+x1(k-2)+0.5*x1(k)^2+2*x1(k)*x1(k-2)-1.6*x1(k-1)^2+1.2*x1(k-2)^2+0.8*x1(k-1)*x1(k-2);

end

Pn = 0.01;
d = d+sqrt(Pn)*n;


%%identification
% h1 = [1, -0.7, -0.5, 0.4];
% x  = rand(1100,1);
% x_aux =  flipud(buffer(x,4,4-1));
% X1 = h1*x_aux;
% f  = @(x) tanh(x);
% Y  = f(X1);
% Pn = 10^(-20/10)*var(Y);
% d  = Y + sqrt(Pn)*randn(1,1100); 

%% prediction
% u1 = uu((ee-1)*1100+1:ee*1100);
% x = u1;
% d = uu((ee-1)*1100+1+L:ee*1100+L);

% input = data_x((ee-1)*1098+1:ee*1098);
% x = input + 0.0001;
% d = x(1+L:end); 
% x = x(1:end-L);

%% equalization
% u = randn(1,1100+N-1) > 0;
% u = 2*u - 1;
% d = double(u(L+1:L+1100)); 
% % nonlinear channel
% y = u + 0.5*[0, u(1:end-1)];        % linear channel
% y = y + 0.2*y.^2+ 0.1*y.*3;                   % nonlinear channel
% x = y + sqrt(0.001)*randn(1,length(u));     % channel noise

x_aux          =   x;
x_aux          =   buffer(x_aux,N,N-1,'nodelay');%equalization 'nodelay'
x_tdl          =   flipud(x_aux);
    
[~,M]          =   size(x_tdl);

prob = 0.0;
bernp1 = rand(1,1000) <= prob;

oi = sqrt(25)*randn(1,K);
d(bernp1==1)=oi(bernp1==1);

cond1=0;
cond2=0;

P_up_all = 0.001:0.1*(1-0.001):1;

g=zeros(N,1);
c=zeros(N,1);

for p = 1:length(P_up_all)
    
P_up =  P_up_all(p);                
beta = ones(1,1991);

X       = x_tdl(:,1);
q{1}  = km_kernel(X,X,'gauss',sigma_k);
G_t = q{1};
G_M = G_t;
eta{1}  = d(1)/q{1};
e{1} = 0;
y(1) = d(1);
G1      = 1;
d_til   = d(1);
x_tdl_til = x_tdl(:,1);

threshold_tau = floor(0.2*size(x_aux,2));
teste1 = [];
    e_antigo = 0;
b = 0.99;
c = 1;

for k=1:1000

% tau_max = 40000;
%   if  k < threshold_tau
%       sqrt_tau_max1 = inf;
%       sqrt_tau_max2 = -inf;
%   end
%   if k >= threshold_tau
%       sqrt_tau_max1 = mean(teste1) + 3*sqrt(var(teste1));
%       sqrt_tau_max2 = mean(teste1) - 3*sqrt(var(teste1));
%   end  
sqrt_tau_max1 = sqrt(40000000000000000000000);
sqrt_tau_max2 = -sqrt(40000000000000000000000);

   alpha2 = P_up*N*(1-lambda)/(2-P_up*(1-lambda));

   sqrt_tau = sqrt(1+alpha2)*qfuncinv((P_up)/2);
%    sqrt_tau = sqrt(1+alpha2)*qfuncinv((P_up+2*qfunc(sqrt_tau_max1/sqrt(1+alpha2)))/2);

    eee(k) = d(k+1)' - km_kernel(X,x_tdl(:,k+1),'gauss',sigma_k)'*eta{c};
    y(k+1) = km_kernel(X,x_tdl(:,k+1),'gauss',sigma_k)'*eta{c};
    e_atual = abs(eee(k)^2);
    teste(k) = abs(eee(k))/sqrt((1-b)*e_atual+b*e_antigo);
%     teste2(k) = (eee(k))/sqrt(0.001);
    
    if (teste(k) <= sqrt_tau)
        delta(k)=0;
        cond1=cond1+1;
        
    elseif (teste(k) > sqrt_tau_max1) || teste(k) < sqrt_tau_max2
        delta(k)=0;
        cond2 = cond2+1;
        
    else
        delta(k)=1;
    end
   
    if (delta(k) == 0)
        alpha(k) = 0;
        if (teste(k) > sqrt_tau_max1) || teste(k) < sqrt_tau_max2  
            eee(k) = 0.0;
%             d(k) = 0;
        end
    else
    
    teste1 = [teste1 eee(k)/sqrt(0.01)];    
          
    clear r
    clear v
    d_til = [d_til(1:c) d(k+1)];
    x_tdl_til = [x_tdl_til(:,1:c) x_tdl(:,k+1)];
    e{k}  = d_til' - km_kernel(X,x_tdl_til,'gauss',sigma_k)'*eta{c}; 
    g_k  = km_kernel(X,x_tdl(:,k+1),'gauss',sigma_k);
    X = [X x_tdl(:,k+1)];
    q_k  = km_kernel(x_tdl(:,k+1),x_tdl(:,k+1),'gauss',sigma_k);
    G_t     = [G_t g_k ; g_k' q_k];
    G1    = [lambda*G1;1];
    G_M   = (G1*G1').*G_t;
    r(:,1) = e{k};%%%%%%%%%%5
    gamma(1) = dot(G_M*r(:,1),r(:,1));
    v(:,1)   = G_M*r(:,1);
    beta(1) = 1;
    u       = zeros(2,1);
    
    for j = 1:2
        alpha(j) = eta1*gamma(j)/dot(v(:,j),v(:,j));
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
    eta{c+1} = [eta{c}; 0] + ((G1)).*((r(:,1:end-1))*u);
    c        = c + 1;
    end
     e_antigo = (1-b)*e_atual+b*e_antigo;

end
MSE(:,p)= MSE(:,p)+(abs(eee.')).^2;
P_up_est(p) = size(delta(delta ==1),2)/(size(delta,2));
% misalig_total(:,p) = misalig_total(:,p)+misalig.';
% e_all(:,p) = e;
% tau(p) = sqrt_tau;
d_all(:,p) = d;
% dd_all(:,p) = di;
y_all(:,p) = y;
% if p == 6
%     w_all = w;
% end

end
% w_alll = w_alll + w_all;
% P_uppp(:,ee)=P_up_est;
ee
end

% P_uppp = mean(P_uppp');
MSE = MSE/ensemble;
% misalig_total = misalig_total/ensemble;

figure,
plot(10*log10(MSE));
set(gca,'fontsize',18)
xl = xlabel('Number of iterations, $k$'); 
yl = ylabel('Misalignment (dB)'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% saveas(gcf,'id_data_set2','fig');


figure,
plot(P_up_est,'-ro')
hold on
plot(P_up_all,'-b+')
set(gca,'fontsize',18)
xl = xlabel('Trial number'); yl = ylabel('$\hat{P}_{up},P_{up}$'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
leg1 = legend('$\hat{P}_{up}$','$P_{up}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18); 
%saveas(gcf,'N_1_P_up','fig');

% Plotando as amostras originais e estimadas/previstas:
figure
plot(d_all(:,11),'-bo')
hold on
plot(y_all(:,11),'-+r')
set(gca,'fontsize',18)
xlim([0 1000])
xl = xlabel('Number of iterations, $k$'); yl = ylabel('Signals'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
leg1 = legend('${d(k)}$','$y(k)$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18);
saveas(gcf,'x_track_ds','fig');