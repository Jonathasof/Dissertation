clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Problema Mackey-Glass  - KCG                   %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('data.mat')
data_x = data_wind;
% % data_x = (data_x - mean(data_x))/var(reshape(data_x,65536*5,1))^(1/2);

N = 10;
K       = 5000;  
sigma_k = sqrt(0.1);
L = 3;
SNR = 20;

u = Mackey_glass;
u = u + sqrt(0.001)*randn(size(u));

for ee=1:1
   
%% testes
% x = sqrt(0.001)*randn(1,1100);
% n = sqrt(0.001)*randn(1,1100);
% 
% d(1) = -0.76*x(1) + 0.5*x(1)^2;
% d(2) = -0.76*x(2) - 1*x(1) + 0.5*x(2)^2 - 1.6*x(1)^2;
% for k = 3:1100-L
% 
%     d(k) = -0.76*x(k)-x(k-1)+x(k-2)+0.5*x(k)^2+2*x(k)*x(k-2)-1.6*x(k-1)^2+1.2*x(k-2)^2+0.8*x(k-1)*x(k-2)+n(k);
% 
% end
% num = 1;
% den = [1 -0.95];
% x   = filter(num,den,x);


%% Prediction
u1 = u((ee-1)*1100+1:ee*1100);
x = u1;
d = u((ee-1)*1100+1+L:ee*1100+L);

% input = data_x((ee-1)*1098+1:ee*1098);  %%falta detalhes
% x = input+0.001;
% d = x(1+L:end); 
% x = x(1:end-L);

%% identification

% nn = sqrt(1-0.88^2)*randn(1,K);
% x_0(1) = nn(1);
% for i=2:K
% x_0(i) = 0.88*x_0(i-1)+nn(i);
% end
% n_n =  sqrt(.001)*randn(1,K);
% % h1 = [1, -0.7, -0.5, 0.4];
% h1 = [0.1010 0.3030 0 -0.2020 -0.4040 -0.7071 -0.4040 -0.2020];
% %x  = rand(1100,1);
% 
% x         =   [zeros(N-1,1);(x_0')];
% x_aux=buffer(x,N,N-1,'nodelay');
% x_k=flipud(x_aux);
% 
% di = h1*x_k+n_n;
% d=di;
% f  = @(x) (x);
%Pn = 10^(-20/10)*var(Y);

%% equalization  %%%%falta detalhes  
% u = randn(1,1100+N-1) > 0;
% u = 2*u - 1;
% d = double(u(L+1:L+1100)); 
% % nonlinear channel
% y = u + 0.5*[0, u(1:end-1)];        % linear channel
% y = y + 0.4*y.^2 + 0.05*y.*3;                   % nonlinear channel
% x = y + sqrt(0.001)*randn(1,length(u));     % channel noise

S   =   struct('filterOrderNo',N,'iterationsnumbers',K,'sigma_k',sigma_k);
            
[e] = KCG_exponential_2(x,d,S);

% for i=1:990
% MSE(i,ee) = abs(e{i}(i)).^2;
% end

MSE(:,ee) = abs(e');

end


p = mean(MSE,2);
figure(1)
plot(10*log10(p))
xl = xlabel('$k$'); yl = ylabel('MSE (dB)'); 
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);