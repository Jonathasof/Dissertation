%    clear all; clc

images = loadMNISTImages('train-images.idx3-ubyte');
labels = loadMNISTLabels('train-labels.idx1-ubyte');

images_train = images(:,1:60000);
labels_train = labels(1:60000,:);


 % Test set
images_test = loadMNISTImages('t10k-images-idx3-ubyte');
labels_test = loadMNISTLabels('t10k-labels-idx1-ubyte');


% load('matlab/emnist-letters')
% images_train = dataset.train.images'/255;
% labels_train = dataset.train.labels;
% 
% images_test = dataset.test.images'/255;
% labels_test = dataset.test.labels;


lay = 10;
y_one_hot_train = zeros(size(labels_train,1),lay);
for i = 1:lay
    rows = labels_train == i-1;
    y_one_hot_train( rows, i ) = 1;
end
y_one_hot_train= y_one_hot_train';

y_one_hot_test = zeros(size(labels_test,1),lay);
for i = 1:lay
    rows = labels_test == i-1;
    y_one_hot_test(rows, i) = 1;
end
y_one_hot_test = y_one_hot_test';

K = 128;
learning_rate = 0.1;
layers = [784 1024 1024 1024 lay];%[784 1024 1024 10];
activation_fun = @ReLu;%
der_activation_fun = @der_ReLu;%
output_fun =@Softmax;%
der_output_fun = @Softmax;%
num = 1;
epoch = 100;
loss_fun = 1;

[W, loss_train, loss_test] = Neural_Network_dropout_clas(images_train,y_one_hot_train,layers,learning_rate,activation_fun,...
                der_activation_fun,output_fun,der_output_fun,num,loss_fun,epoch, K, labels_train, labels_test, images_test, y_one_hot_test);
            
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

