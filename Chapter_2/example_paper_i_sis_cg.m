clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         System identification (DS-CG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ensemble = 200;               % number of independent runs
N     = 8;                   % number of coefficients
K=5000;                      % number of iterations
gamma = 10^(-4);                  % regularization factor



lambda = 0.98;
eta = 0.48;
delta2 = 10^(-4);
% sigma_ne = sqrt(1-0.88^2);
MSE = zeros(K,11);
misalig_total = zeros(K,11);
w_alll = zeros(8,5001);

for ee=1:ensemble

% nn = sqrt(1-1.22^2)*randn(1,K);
nn = sqrt(0.25)*randn(1,K);
x_0(1) = nn(1);
x_0(2) = -0.55*x_0(1)+nn(2);
x_0(3) = -0.55*x_0(2)-1.221*x_0(1)+nn(3);
x_0(4) = -0.55*x_0(3)-1.221*x_0(2)-0.4955*x_0(1)+nn(4);
for i=5:K
% x_0(i) = 0.88*x_0(i-1)+nn(i);
x_0(i) = -0.55*x_0(i-1)-1.221*x_0(i-2)-0.49955*x_0(i-3)-0.4536*x_0(i-4) + nn(i);
end
n_n =  sqrt(.001)*randn(1,K); 

h = [0.1010 0.3030 0 -0.2020 -0.4040 -0.7071 -0.4040 -0.2020];


prob = 0.01;
bernp1 = rand(1,K) <= prob;

% System:
%  num=zeros(1,N+1);
%  num(1)=1;
%  num(N+1)=-1;
%  
%  den=[1 -1];
%  d=filter(num,den,x_0)+n_n; 
 %d = [0 d];

x_aux          =   [zeros(N-1,1);(x_0')];
x_aux=buffer(x_aux,N,N-1,'nodelay');
x_k=flipud(x_aux);

di = h*x_k+n_n;
d=di;

oi = sqrt(5)*randn(1,K);
d(bernp1==1)=oi(bernp1==1);
w(:,size(x_aux,2))=zeros(N,1);

cond1=0;
cond2=0;

P_up_all = 0.3;%0.001:0.1*(1-0.001):1;
sigma_n= ones(size(x_aux,2),1);

R=eye(N);

g=zeros(N,1);
c=zeros(N,1);

for p = 1:length(P_up_all)
    
P_up =  P_up_all(p);                

%alpha2 = ((L_r+1)*P_up*mu)/((2-P_up*mu)*(1+L_r*(1-P_up*mu)^2));


threshold_tau = floor(0.02*size(x_aux,2));
teste1 = [];

for k=1:size(x_aux,2)


% tau_max = 40000;
  if  k < threshold_tau
      sqrt_tau_max1 = inf;
      sqrt_tau_max2 = -inf;
  end
  if k >= threshold_tau
      sqrt_tau_max1 = mean(teste1) + 3*sqrt(var(teste1));
      sqrt_tau_max2 = mean(teste1) - 3*sqrt(var(teste1));
  end  
% sqrt_tau_max1 = sqrt(40000000000000000000000);
% sqrt_tau_max2 = -sqrt(40000000000000000000000);

   alpha2 = P_up*N*(1-lambda)/(2-P_up*(1-lambda));

%  sqrt_tau = sqrt(1+alpha2)*qfuncinv((P_up)/2);
   sqrt_tau = sqrt(1+alpha2)*qfuncinv((P_up+2*qfunc(sqrt_tau_max1/sqrt(1+alpha2)))/2);

    y(k) = w(:,k)'*x_k(:,k);
    e(k) = d(k)- y(k);
    
    teste(k) = abs(e(k))/sqrt(0.001);
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
        w(:,k+1)=w(:,k);
%         teste1 = [teste1 e(k)/sqrt(0.001)];    
        if (teste2(k) > sqrt_tau_max1) || teste2(k) < sqrt_tau_max2  
            e(k) = 0.0;
            d(k) = 0;
        end
    else
    
    teste1 = [teste1 e(k)/sqrt(0.001)];    
    
    R = lambda*R+x_k(:,k)*x_k(:,k)';
    
    alpha(k) = eta * (c' * g)/(c' * R * c + delta2);
    
    g1 = g;
    
    g= lambda * g1 - alpha(k) *  R*c + x_k(:,k) * e(k);
    
    w(:,k+1)= w(:,k) + delta(k)*alpha(k) * c;
    
    beta(k) = ((g - g1)' * g)/(g1' * g1+delta2);
    
    c=g + beta(k) * c;
        
    end
    misalig(:,k) = (norm(w(:,k+1) - h.'))^2/norm(h)^2;
end
MSE(:,p)= MSE(:,p)+(abs(e.')).^2;
misalig_total(:,p) = misalig_total(:,p)+misalig.';
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
% w_alll = w_alll + w_all;
P_uppp(:,ee)=P_up_est;
ee
end

P_uppp = mean(P_uppp');
MSE = MSE/ensemble;
misalig_total = misalig_total/ensemble;

figure,
plot(10*log10(misalig_total(:,1)));
set(gca,'fontsize',18)
xl = xlabel('Number of iterations, $k$'); 
yl = ylabel('Misalignment (dB)'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% saveas(gcf,'id_data_set2','fig');


% figure,
% plot(x_0(1:5:end));
% set(gca,'fontsize',18)
% xl = xlabel('Number of iterations, $k$'); 
% yl = ylabel('Samples'); 
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% saveas(gcf,'id_data_set2','fig');
% 
% 
% figure,
% plot(P_uppp,'-ro')
% hold on
% plot(P_up_all,'-b+')
% set(gca,'fontsize',18)
% xl = xlabel('Trial number'); yl = ylabel('$\hat{P}_{\rm up},P_{\rm up}$'); 
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% leg1 = legend('$\hat{P}_{\rm up}$','$P_{\rm up}$');
% set(leg1,'Interpreter','latex');
% set(leg1,'FontSize',18); 
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
% 
% % figure,
% % plot(teste)
% % set(gca,'fontsize',18)
% % xl = xlabel('Number of iterations, $k$'); yl = ylabel('$|e|/\sigma_x$'); 
% % set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% % set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% % saveas(gcf,'N_1_e_bar_sigma','fig');
% 
% % figure,
% % plot(10*log10(MSE(:,6)));
% % xl = xlabel('$k$'); yl = ylabel('MSE (dB)'); 
% % set(xl,'Interpreter','latex');set(xpup=[2 4 5 11];
% 
% 
% names = {'id_x_track_1','id_x_track_2','id_x_track_3','id_x_track_4'};
% 
% % set(yl,'Interpreter','latex');set(yl,'FontSize',14);
% % saveas(gcf,'mse_ds','fig');
% 
% 
% figure
% plot(dd_all(:,6),'-bo')
% hold on
% plot(y_all(:,6),'-+r')
% set(gca,'fontsize',18)
% xl = xlabel('Number of iterations, $k$'); 
% yl = ylabel('$signals$'); 
% set(xl,'Interpreter','latex');set(xl,'FontSize',14);
% set(yl,'Interpreter','latex');set(yl,'FontSize',14);
% leg1 = legend('${d(k)}$','$y(k)$');
% set(leg1,'Interpreter','latex');
% set(leg1,'FontSize',14);
% % saveas(gcf,'x_track_ds','fig');
% 
% 
% pup=[2 4 5 11];
% names = {'id_x_track_1','id_x_track_2','id_x_track_3','id_x_track_4'};
% for j=1:length(pup)
% pp=pup(j);    
% figure
% % plot(d_all(:,pp),'-bo')
% % hold on
% plot(y_all(:,pp),'-+r')
% hold on
% plot(dd_all(:,pp),'m--x')
% set(gca,'fontsize',18)
% % ylim([-3 3])
% xlim([4900 5000])
% xl = xlabel('Number of iterations, $k$'); yl = ylabel('Output signal'); 
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% leg1 = legend('$y(k)$','${d(k)}$');
% set(leg1,'Interpreter','latex');
% set(leg1,'FontSize',18);
% saveas(gcf,names{j},'fig');
% end
% 
% 
% w_ok = (w_alll/ensemble);
% figure,
% for i=1:N
% plot(w_ok(i,:));
% hold on
% end
% set(gca,'fontsize',18)
% ylim([-1 0.6])
% xlim([0 5000])
% line(xlim, [1 1],'Color','black','LineStyle','--');
% xl = xlabel('Number of iterations, $k$'); 
% yl = ylabel('$w_i$, para $i = 1 \ldots 8$','interpreter','latex');
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% saveas(gcf,'id_coef_evol_ds','fig');
% 
% figure,
% plot(10*log10(MSE(:,6)),'r')
% % hold on
% % plot(10*log10(MSE(:,11)))
% xl = xlabel('Number of iterations, $k$'); 
% yl = ylabel('MSE (dB)');
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% % legend('50%')
