clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Predição com seleção de dados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ensemble=1;
N     = 20;                   % número de coeficientes
L = 1;
K=5000;  % número de iterações
gamma = 10^(-4);                  % constante para regularização


sigma_nn = sqrt(0.001);
lambda = 0.98;
eta = 0.48;
delta2 = 10^(-4);

MSE = zeros(K-N+1,11);
for ee=1:ensemble
x = randn(1,K);
x = sign(x);

n = (sigma_nn *randn(1,K));

num = [1 1.6561 0.34];
den = 1;%[1 1 0.81];

x_aux = filter(num,den,x)+n;
x_aux=buffer(x_aux,N,N-1,'nodelay');
x_k=flipud(x_aux);

di = x';
d=di;

w(:,size(x_aux,2))=zeros(N,1);


cond1=0;
cond2=0;

P_up_all = 0.001:0.1*(1-0.001):1;
% P_up_all =0.5;

R=eye(N);

g=zeros(N,1);
c=zeros(N,1);

for p = 1:length(P_up_all)




    
P_up =  P_up_all(p);                % probabilidade de atualização

%alpha2 = ((L_r+1)*P_up*mu)/((2-P_up*mu)*(1+L_r*(1-P_up*mu)^2));


threshold_tau = floor(0.2*size(x_aux,2));
e_antigo = 0;
b = 0.99;

for k=1:size(x_aux,2)

tau_max = 40;
  

%    alpha2 = P_up*N*(1-lambda)/(2-P_up*(1-lambda));

 sqrt_tau = qfuncinv(P_up/2);
  
    y(k) = w(:,k)'*x_k(:,k);
    e(k) = d(k)- y(k);
    e_atual = e(k)^2;
    
    teste(k) = abs(e(k))/sqrt((1-b)*e_atual+b*e_antigo);
    ooi(p,k) = ((1-b)*e_atual+b*e_antigo);
%     teste(k) = abs(e(k))/sigma_nn;
    
    
   
    
    if (teste(k) <= sqrt_tau)
        delta(k)=0;
        cond1=cond1+1;
        
    elseif (teste(k) > (tau_max))
        delta(k)=0;
        cond2 = cond2+1;
        
    else
        delta(k)=1;
    end
   
    if (delta(k) == 0)
        w(:,k+1)=w(:,k);
        if teste > sqrt(tau_max)
            
            e(k) = 0;
            d(k) = 0;
        end
    else
    
    
    
    
    
    R = lambda*R+x_k(:,k)*x_k(:,k)';
   
    
    
    alpha(k) = eta * (c' * g)/(c' * R * c + delta2);
    
    g1 = g;
    
    g= lambda * g1 - alpha(k) *  R*c + x_k(:,k) * e(k);
    
    w(:,k+1)= w(:,k) + delta(k)*alpha(k) * c;
    
    beta(k) = ((g - g1)' * g)/(g1' * g1+delta2);
    
    c=g + beta(k) * c;
        
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
w_all{p} = w;
end


ee
end


MSE = MSE/ensemble;
% 
% figure,
% plot(x_0(1:5:end));
% set(gca,'fontsize',18)
% xl = xlabel('Number of iterations, $k$'); 
% yl = ylabel('Samples'); 
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% saveas(gcf,'id_data_set2','fig');
% 
% % 
figure,
plot(P_up_est,'-ro')
hold on
plot(P_up_all,'-b+')
set(gca,'fontsize',18)
xl = xlabel('Trial number'); yl = ylabel('$\hat{P}_{\rm up},P_{\rm up}$'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
leg1 = legend('$\hat{P}_{\rm up}$','$P_{\rm up}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18); 
% saveas(gcf,'id_P_up','fig');
% 
% 
% figure,
% plot(abs(e_all(:,11)),'-m*')
% hold on
% plot(sigma_n,'-b+')
% set(gca,'fontsize',18)
% xl = xlabel('Number of iterations, $k$'); yl = ylabel('$|e|,\sigma_x$'); 
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% leg1 = legend('$|e|$','$\sigma_x$');
% set(leg1,'Interpreter','latex');
% set(leg1,'FontSize',18);
% saveas(gcf,'N_1_sigma_x_and_e','fig');



d_all(1,:)=[];

d_all(end,:)=[];

pup=[1 2 11];
for j=1:length(pup)
pp=pup(j);    
figure
Q=plot(1:K-N+L,(d_all(1:K-N+L,pp)),'-+m',1:K-N+L,(y_all(1:K-N+L,pp)),'x--r');
set(Q,{'LineWidth'},{2})

set(gca,'fontsize',18)
ylim([-1.5 1.5])
xlim([4950 5000])
xl = xlabel('Number of iterations, $k$'); yl = ylabel('Output signal'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
leg1 = legend('${y(k)}$','$\hat{y}(k)$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18);
end


figure,
plot(10*log10(MSE(:,2)))
hold on
plot(10*log10(MSE(:,11)))
xl = xlabel('Number of iterations, $k$'); yl = ylabel('MSE (dB)'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
leg1 = legend('DS-CG 1', 'DSCG 11');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',12);

figure,
for i=1:11
plot(ooi(i,:))
hold on
end