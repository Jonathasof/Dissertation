function [W, loss_train, acc_train, loss_test, acc_test] = Neural_Network_dropout2(x,y,layers,learning_rate,...
                                        activation_fun,der_activation_fun,...
                                        output_fun,der_output_fun,num,loss_fun,epoch,K, labels_train, labels_test, x_test, y_test, P_up)

% Dropout somente zerando as entradas das matrizes

c = 0;
eee = [];
% The number of training vectors.
trainingSetSize = size(x, 2);
%P_up = 0.5; 
x = [ones(1,size(x,2));x];

% Initialize the weights for the hidden layer and the output layer.  
for i=1:size(layers,2)-1
    W{i} = randn(layers(i+1),layers(i));
    delta_w{i} = zeros(layers(i+1),layers(i)+1);
    delta_w_antigo{i} = zeros(layers(i+1),layers(i)+1);
    W{i} = [0.01*ones(layers(i+1),1) W{i}]; %come�ando com bias zerado
    W{i} = W{i}./sqrt(layers(i)/2);
end 


iter = floor(trainingSetSize/K);
dd = 1;
ee = 1;
% b = 0.9;
% mu_antigo = zeros(10,1);
% sigma_antigo = zeros(10,1);

for ep = 1: epoch
%     H = waitbar(0, 'Please wait...');
    c = 0;
    P_up_est(ep) = 0;
%     sigma_antigo = P_up*sigma_antigo;
%     mu_antigo = P_up*mu_antigo;
    
    for t = 1: iter
%         if ~ishghandle(H)
%            error('User closed waitbar.');
%         end
        c = c + 1;
%         H = waitbar(c / iter, H, sprintf('%d of %d', c, iter));

         % Select which input vector to train on.
         n = randperm(size(x,2),K);
            
         % Propagate the input vector through the network.
         inputVector = x(:, n);
         s{1} = double(inputVector);
         
         for r=1:size(layers,2)-2 
            e{r} = W{r}*s{r};
            s{r+1} = [ones(1,K); activation_fun(e{r})];
         end
         
         outputActualInput = W{size(layers,2)-1}*s{size(layers,2)-1};
         
         for nn=1:K
              outputVector(:,nn) = output_fun(outputActualInput(:,nn));  
              der_outputActualInput(:,nn) = der_output_fun(outputActualInput(:,nn)); 
         end
         
         targetVector = y(:,n);
            
         % Backpropagate the errors.
         if num == 1
             delta{size(layers,2)-1} = (outputVector - targetVector); %cross entropy com softmax
         elseif num == 2
             delta{size(layers,2)-1} = der_outputActualInput.*(outputVector - targetVector);
         elseif num == 3
             delta{size(layers,2)-1} = (outputVector.*targetVector); % relative entropy com sigmoid
         else
             print('Fun��o n�o existe')
             break
         end
         

%           eee{3} = -targetVector.*log(outputVector)-(1-targetVector).*log(1-outputVector);
          eee{3} = (outputVector - targetVector).^2;

%          cc2 = [cc2 sum(eee{3})];
%          histogram(cc2/10,'Normalization','pdf')     

         nn = binornd(K,P_up);
         if nn == 0
             nn = 1;
         end
         nn_all = sum(eee{3});
         nn_all_sort = sort(nn_all);
         n22 = nn_all >= nn_all_sort(K+1 - nn); %>= finv(P_up,3.25,300);
         P_up_est(ep) = P_up_est(ep) + sum(n22)/K;
            
         delta{end} = delta{end}(:,n22);

         for r=size(layers,2)-2:-1:1
             er{r} = W{r+1}(:,2:end)'*delta{r+1};
             delta{r} = der_activation_fun(e{r}(:,n22)).*(er{r});  
         end 
 
         for i = 1:size(layers,2)-2
             delta_w{i+1} =  - (learning_rate/(sum(n22)))*delta{i+1}*s{i+1}(:,n22)';
             W{i+1} = W{i+1} + delta_w{i+1};
         end
         delta_w{1} =  - (learning_rate/(sum(n22)))*delta{1}*s{1}(:,n22)';
         W{1} = W{1} + delta_w{1};

    end
%     b = b*0.99;
      P_up_est(ep) = P_up_est(ep)/iter; 
         
        % Calculate the loss and accuracy for plotting.
         s{1} = double(x);
         s_test{1} = double(x_test);
         s_test{1} = [ones(1,size(s_test{1},2));s_test{1}];
         
         targetVector = y; 
         targetVector_test = y_test;
         for r=1:size(layers,2)-2 
             if r == 1
                 e{r} = 1*W{r}*s{r};
                 z_test{r} = 1*W{r}*s_test{r};
             else
                 e{r} = dd*W{r}*s{r};
                 z_test{r} = dd*W{r}*s_test{r};
             end
             s{r+1} = [ones(1,size(x,2)); activation_fun(e{r})];
             s_test{r+1} = [ones(1,size(x_test,2)); activation_fun(z_test{r})];
         end
         outputActualInput = dd*(W{size(layers,2)-1})*s{size(layers,2)-1};
         for n=1:size(x,2)
              outputVector_train(:,n) = output_fun(outputActualInput(:,n));   
         end
%          outputVector = output_fun(outputActualInput);
         outputActualInput_test = dd*(W{size(layers,2)-1})*s_test{size(layers,2)-1};
         for n=1:size(x_test,2)
              outputVector_test(:,n) = output_fun(outputActualInput_test(:,n));   
         end
%          outputVector_test = output_fun(outputActualInput_test);
         [~,i] = max(outputVector_train',[],2);
         acc_train(ep) = sum((i-1) == labels_train)/size(x,2);
         [~,i] = max(outputVector_test',[],2);
         acc_test(ep) = sum((i-1) == labels_test)/size(x_test,2);
             
         if loss_fun == 1 %cross entropy
             loss_train(ep) =  - 1/size(x,2)*sum(sum(targetVector.*log(outputVector_train+0.00001)+(1-targetVector).*log(1-outputVector_train+0.00001)));
             loss_test(ep) =  - 1/size(x_test,2)*sum(sum(targetVector_test.*log(outputVector_test+0.00001)+(1-targetVector_test).*log(1-outputVector_test+0.00001)));
         elseif loss_fun == 2 % squared error
             loss_train(ep) = 1/(2*size(x,2))*sum((vecnorm((outputVector_train - targetVector)', 2)));
             loss_test(ep) = 1/(2*size(x_test,2))*sum((vecnorm((outputVector_test - targetVector_test)', 2)));
         elseif loss_fun == 3 % Relative entropy
             loss_train(ep) = - 1/size(x,2)*sum(sum(targetVector.*log(outputVector_train./(targetVector+0.0001))));
             loss_test(ep) = - 1/size(x_test,2)*sum(sum(targetVector_test.*log(outputVector_test./(targetVector_test+0.0001))));             
         end
% learning_rate = learning_rate*1/(1+1e-6*iter); 
% close(H)
fprintf('%d/%d %.4f %.4f %.4f %.4f %.6f\n',...
         ep, epoch, loss_train(ep), acc_train(ep), loss_test(ep), acc_test(ep), P_up_est(ep));
     
end


end