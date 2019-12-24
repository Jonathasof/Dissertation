clear all;
close all;

erro = [];
w = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Predição com seleção de dados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('data.mat')

% filename = 'AirQualityUCI.xlsx';
% A = xlsread(filename);

N = 8;                   % número de coeficientes
L = 1;                   % numero de amostras a frente
%N_samples=size(A,1)-1;  % número de iterações
L_r = 2;                 % fator de reuso de dados

lambda = 0.98;
eta = 0.48;
delta2 = 10^(-7);
data_x = data_wind;
% data_x = A(2:end,13);
% data_x(data_x==-200)=[];
% N_samples=size(data_x,1)-1;
% input = data_x(1:N_samples); 
% x = input;


MSE = zeros(8192-N-L+1,11);
cc = 0;
ensemble = 40;

for mc = 1:ensemble
    
input = data_x((mc-1)*8192+1:mc*8192);
cc = cc + 1

x = input + 0.0001;
x_aux=buffer(x,N,N-1,'nodelay');
x_k=flipud(x_aux);
w(:,size(x_aux,2)-L_r)=zeros(N,1);
d = x_k(1,L+1:end);

cond1=0;
cond2=0;

P_up_all = 0.001:0.1*(1-0.001):1;
%P_up_all = 0.0001:0.1:1.0001;
sigma_n= ones(size(x_aux,2)-L,1);

R=eye(N);

g=zeros(N,1);
c=zeros(N,1);

for p = 1:length(P_up_all)

P_up =  P_up_all(p); % probabilidade de atualização


%alpha1 = ((L_r+1)*P_up*mu)/((2-P_up*mu)*(1+L_r*(1-P_up*mu)^2));%AP
% alpha1 = P_up*(N + 1)*(1-lambda)/(2-P_up*(1-lambda)); %RLS
%alpha1= lambda*eta*P_up*mu/(2-lambda*eta*P_up*mu);%Proposta

threshold_tau = floor(0.2*size(x_aux,2)-L);


for k=1:size(x_aux,2)-L

  if  k < threshold_tau
      tau_max = inf;
  end
  if k >= threshold_tau
      tau_max = mean(teste) + 3*sqrt(var(teste));
  end
  
  r(1) = x(k)*x(k);
  beta = 0.01;

 if k>N+1
for l=2:N
    r(l)=(beta)*x(k)*x(k-l)+(1-beta)*r(l-1);
end
xi_min = r(1)-w(:,k)'*r';

sigma_n(k) = sqrt(abs(xi_min));
 end
 
    sqrt_tau = qfuncinv(P_up/2);
%     sqrt_tau = qfuncinv((P_up+2*qfunc(tau_max))/2);

    y(k) = (x_k(:,k)'*w(:,k));
    e(k)=d(k)-y(k);
    
    teste(k) = abs(e(k))/sigma_n(k);

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
        if teste > (tau_max)
            e(k) = 0;
            d(k) = 0;
        end
    else
       
    R = lambda*R+x_k(:,k)*x_k(:,k)';
    
    alpha = eta * (c' * g)/(c' * R * c + delta2);
    
    g1 = g;
    
    g= lambda * g1 - alpha *  R*c + x_k(:,k) * e(k);
    
    w(:,k+1)= w(:,k) + alpha * c;
    
    beta = ((g - g1)' * g)/(g1' * g1+delta2);
    
    c = g + beta * c;
    
    end 
        
    
    
end
MSE(:,p)= MSE(:,p)+(abs(e.')).^2;
P_up_est(p) = size(delta(delta ==1),2)/(size(delta,2));
e_all(:,p) = e;
tau(p) = sqrt_tau;
d_all(:,p) = d;
y_all(:,p) = y;
w_all{p} = w;
end
P_uppp(:,mc)=P_up_est;

end
P_uppp = mean(P_uppp');
% for i=1:length(P_up_all)
%     erro_med(i) = mean ((abs(e_all(8500:8950,i))).^2);
% end

% erro = [erro; erro_med];
% 
%  end 
  
% oi=data_x;
% data_x = oi(8:1:31);

% figure,
% plot(0:23,data_x);
% set(gca,'fontsize',18)
% xlim([0 23])
% xl = xlabel('Time $(h)$'); yl = ylabel('Temperature in $^{\circ}$C'); 
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% %saveas(gcf,'data_set2','fig');


figure,
plot(P_up_all,'-b+')
hold on
plot(tau,'-mo')
leg1 = legend('$P_{up}$','$\sqrt{\tau}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18);
%saveas(gcf,'tau','fig');

figure,
plot(P_uppp,'-ro')
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


figure,
plot(abs(e_all(:,11)),'-m*')
hold on
plot(sigma_n,'-b+')
set(gca,'fontsize',18)
xlim([0 9000])
xl = xlabel('Number of iterations, $k$'); yl = ylabel('$|e|,\sigma_x$'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
leg1 = legend('$|e|$','$\sigma_x$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18);
%saveas(gcf,'N_1_sigma_x_and_e','fig');

figure,
plot(teste)
set(gca,'fontsize',18)
xlim([0 9000])
xl = xlabel('Number of iterations, $k$'); yl = ylabel('$|e|/\sigma_x$'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
%saveas(gcf,'N_1_e_bar_sigma','fig');

MSE = MSE/ensemble;

figure,
plot(10*log10(MSE(:,6)));
xlim([0 9000])
xl = xlabel('Number of iterations, $k$'); yl = ylabel('MSE (dB)'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
saveas(gcf,'mse_ds','fig');

% Plotando as amostras originais e estimadas/previstas:
figure
plot(d_all(:,6),'-bo')
hold on
plot(y_all(:,6),'-+r')
set(gca,'fontsize',18)
xlim([0 9000])
xl = xlabel('Number of iterations, $k$'); yl = ylabel('Signals'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
leg1 = legend('${d(k)}$','$y(k)$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18);
%saveas(gcf,'x_track_ds','fig');


pup=[2 5 8 11];
names = {'N_1_x_track_1','N_1_x_track_2','N_1_x_track_3','N_1_x_track_4'};
for j=1:length(pup)
pp=pup(j);    
figure
plot(d_all(:,pp),'-bo')
hold on
plot(y_all(:,pp),'-+r')
set(gca,'fontsize',18)
xlim([8000 8192])
xl = xlabel('Number of iterations, $k$'); yl = ylabel('Predictor signal'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',18);
set(yl,'Interpreter','latex');set(yl,'FontSize',18);
leg1 = legend('${d(k)}$','$y(k)$');
% leg1 = legend('${y(k)}$','$\hat{y}(k)$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18);
%saveas(gcf,names{j},'fig');
end

% w_ok = (w_all{6});
% figure,
% for i=1:N
% plot(w_ok(i,:));
% hold on
% end
% set(gca,'fontsize',18)
% line(xlim, [1 1],'Color','black','LineStyle','--');
% xlim([0 9000])
% xl = xlabel('Number of iterations, $k$'); yl = ylabel('$w_1$','interpreter','latex');
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% %saveas(gcf,'coef_evol_ds','fig');


% figure,
% plot(1:50,erro(:,11),'-*');
% set(gca,'fontsize',18)
% xlim([0 50])
% xl = xlabel('Order, $N$'); yl = ylabel('MSE'); 
% set(xl,'Interpreter','latex');set(xl,'FontSize',18);
% set(yl,'Interpreter','latex');set(yl,'FontSize',18);
% saveas(gcf,'NvsMSE','fig');
