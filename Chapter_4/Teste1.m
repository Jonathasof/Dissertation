% clear all; clc
load('matlab/music')
x = video(:,[1:17,19:end]);
x = (x - min(x))./(max(x) - min(x));
images_train = x(1:55000,:)';
images_test  = x(55001:end,:)';
y = video(:,18);
y = (y - min(y))/(max(y) - min(y));
y_one_hot_train = y(1:55000);
y_one_hot_test  = y(55001:end);
%P_up = 0.3
%% fifa 2018
% load('matlab/X')
% load('matlab/Y')
% Y = (Y-min(Y))/(max(Y)-min(Y));
% images_train = X(1:12000,:)';
% y_one_hot_train = Y(1:12000);
% images_test  = X(12001:end,:)';
% y_one_hot_test  = Y(12001:end);
% labels_train =  Y(1:12000);
% labels_test = Y(12001:end);

%% SuperCond
% load('matlab/X_cond')
% load('matlab/Y_cond')
% el_corr = [19 21 31 49 52 70]; %corr > 0.1 (corr/max(corr))
% shuffle = randperm(length(x));
% train   = x(shuffle,:);
% train   = train(:,[1:18,20,22:30,32:48,50,51,53:69,71:end]); 
% y       = y(shuffle); 
% for i=1:73
%     train(:,9+i) = (train(:,9+i) - min(train(:,9+i)))/(max(train(:,9+i))-min(train(:,9+i)));
% end
% y = (y-min(y))/(max(y)-min(y));
% images_train = train(1:18000,:)';
% y_one_hot_train = y(1:18000);
% images_test  = train(18001:end,:)';
% y_one_hot_test  = y(18001:end);

%% Online Share
% load('matlab/X_online1')
% load('matlab/Y_online1')
% el_corr = [3  23 25:29 55]; %corr > 0.05 
% x = X_online1;
% y = Y_online1;
% shuffle = randperm(length(x));
% x  = x(:,[1:2,4:22,24,30:54,56:end]); 
% x   = x(shuffle,:);
% y  = y(shuffle); 
% x = (x - min(x))./(max(x)-min(x));
% y = (y-min(y))/(max(y)-min(y));
% images_train = x(1:35000,:)';
% y_one_hot_train = y(1:35000);
% images_test  = x(35001:end,:)';
% y_one_hot_test  = y(35001:end);


%% Facebook
% load('matlab/X_fb')
% load('matlab/Y_fb')
% X = X(:,[1:34,36,37,40:end]);
% y = (Y-min(Y))/(max(Y)-min(Y));
% x = (X - min(X))./(max(X)-min(X));
% 
% images_train = x(1:40949,:)';
% y_one_hot_train = y(1:40949);
% images_test  = x(40949+1:end,:)';
% y_one_hot_test  = y(40949+1:end);


K = 2*256;%2*256;
learning_rate = 0.1;
layers = [25 128 128 128 1];%[784 1024 1024 10];
activation_fun = @ReLu;%
der_activation_fun = @der_ReLu;%
output_fun = @Linear;%@Softmax;% 
der_output_fun = @der_Linear;%@der_Softmax;%
num = 2;
epoch = 200;
loss_fun = 2;

 %% Test set
% % images_test = loadMNISTImages('t10k-images-idx3-ubyte');
% % labels_test = loadMNISTLabels('t10k-labels-idx1-ubyte');

%images_test = dataset.test.images'/255;
%labels_test = dataset.test.labels;
% labels_test = labels_test + 1;

% % y_one_hot_test = zeros(size(labels_test,1),10);
% % 
% % for i = 1:10
% %     rows = labels_test == i-1;
% %     y_one_hot_test(rows, i) = 1;
% % end
% % y_one_hot_test = y_one_hot_test';

% [y_hat_test,loss_test] = forward_loss(images_test,y_one_hot_test,10000,layers,W,activation_fun,output_fun,loss_fun);

[W, loss_train, loss_test] = Neural_Network_dropout(images_train,y_one_hot_train,layers,learning_rate,activation_fun,...
                der_activation_fun,output_fun,der_output_fun,num,loss_fun,epoch, K, images_test, y_one_hot_test, 0.8);
            
% [y_hat_train,~] = forward_loss(images_train,y_one_hot_train',60000,layers,W,activation_fun,output_fun,loss_fun);

% figure(1)
% plot(error)
% [~,i] = max(y_hat_train,[],2);
% sum((i-1)== labels_train)/60000
% C  = confusionmat(i-1,labels_train);
% c = {'0';'1';'2';'3';'4';'5';'6';'7';'8';'9'};
% figure(2)
% heatmap(c,c,C);
% figure(3)
% display_network(images(:,1:100)); % Show the first 100 images


% [~,i] = max(y_hat_test,[],2);
% sum((i-1)==labels_test)/10000

