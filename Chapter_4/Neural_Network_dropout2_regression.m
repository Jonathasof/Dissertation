function [W, loss_train, loss_test] = Neural_Network_dropout2_regression(x,y,layers,learning_rate,...
                                        activation_fun,der_activation_fun,...
                                        output_fun,der_output_fun,num,loss_fun,epoch,K, x_test, y_test,P_up)

% The number of training vectors.
trainingSetSize = size(x, 2);
x = [ones(1,size(x,2));x];

% Initialize the weights for the hidden layer and the output layer.  
for i=1:size(layers,2)-1
    W{i} = randn(layers(i+1),layers(i));
    delta_w{i} = zeros(layers(i+1),layers(i)+1);
%     delta_w_antigo{i} = zeros(layers(i+1),layers(i)+1);
    W{i} = [0.0*ones(layers(i+1),1) W{i}]; %come�ando com bias zerado
    W{i} = W{i}./sqrt(layers(i)/2);
end 

iter = floor(trainingSetSize/K);

dd = 1;
b = 0.995;
cc = 1;
sigma_antigo = 0;

for ep = 1: epoch
%     H = waitbar(0, 'Please wait...');
%     c = 0;
    P_up_est(ep) = 0;
    cc2 = [];
    sigma_antigo = P_up*sigma_antigo;

    for t = 1: iter
%         if ~ishghandle(H)
%            error('User closed waitbar.');
%         end
%         c = c + 1;
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
         
         targetVector = y(n,:);
         targetVector = targetVector';
            
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
         
         sigma_atual = (var(outputVector-targetVector));
         eee{3} = (outputVector-targetVector)./(cc*sqrt(((1-b)*sigma_atual+b*sigma_antigo)));%Normalizando
         sigma_antigo = (1-b)*sigma_atual+b*sigma_antigo;
         
         sqrt_tau = qfuncinv((P_up)/2);
         cc2 = [cc2 eee{3}];
         uu = abs(eee{3}) >= sqrt_tau;
         P_up_est(ep) = P_up_est(ep) + sum(uu)/K;

         delta{size(layers,2)-1} = delta{size(layers,2)-1}(:,uu);
         for r=size(layers,2)-2:-1:1
             er{r} = W{r+1}(:,2:end)'*delta{r+1};
             delta{r} = der_activation_fun(e{r}(:,uu)).*(er{r});  
         end 
 
         for i = 1:size(layers,2)-2
             delta_w{i+1} =  - (learning_rate/(sum(uu)))*delta{i+1}*s{i+1}(:,uu)';
             W{i+1} = W{i+1} + delta_w{i+1};
         end
         delta_w{1} =  - (learning_rate/(sum(uu)))*delta{1}*s{1}(:,uu)';
         W{1} = W{1} + delta_w{1};

    end
     P_up_est(ep) = P_up_est(ep)/iter; 
    
        % Calculate the loss and accuracy for plotting.
         s{1} = double(x);
         s_test{1} = double(x_test);
         s_test{1} = [ones(1,size(s_test{1},2));s_test{1}];
         
         targetVector = y'; 
         targetVector_test = y_test';
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
         outputActualInput_test = (W{size(layers,2)-1})*s_test{size(layers,2)-1};
         for n=1:size(x_test,2)
              outputVector_test(:,n) = output_fun(outputActualInput_test(:,n));   
         end
           
         
         r_train = 1 - sum((outputVector_train - targetVector).^2)/sum((targetVector - mean(targetVector)).^2);
         r_test = 1 - sum((outputVector_test - targetVector_test).^2)/sum((targetVector_test - mean(targetVector_test)).^2);
         
         if loss_fun == 1 %cross entropy
             loss_train(ep) =  - 1/size(x,2)*sum(sum(targetVector.*log(outputVector_train+0.00001)+(1-targetVector).*log(1-outputVector_train+0.00001)));
             loss_test(ep) =  - 1/size(x_test,2)*sum(sum(targetVector_test.*log(outputVector_test+0.00001)+(1-targetVector_test).*log(1-outputVector_test+0.00001)));
         elseif loss_fun == 2 % squared error
             loss_train(ep) = 1/(size(x,2))*sum((outputVector_train - targetVector).^2);
             loss_test(ep) = 1/(size(x_test,2))*sum((outputVector_test - targetVector_test).^2);
         elseif loss_fun == 3 % Relative entropy
             loss_train(ep) = - 1/size(x,2)*sum(sum(targetVector.*log(outputVector_train./(targetVector+0.0001))));
             loss_test(ep) = - 1/size(x_test,2)*sum(sum(targetVector_test.*log(outputVector_test./(targetVector_test+0.0001))));             
         end
% learning_rate = learning_rate*1/(1+1e-6*iter); 
% R_squared_train(ep) = 1 - sum((outputVector_train - targetVector).^2)/sum((mean(targetVector) - targetVector).^2);
% R_squared_test(ep) = 1 - sum((outputVector_test - targetVector_test).^2)/sum((mean(targetVector_test) - targetVector_test).^2);
% b = b;
         histogram(cc2,'Normalization','pdf')
         hold on
         plot(-5:0.01:5,normpdf(-5:0.01:5))
         hold off
% close(H)
fprintf('%d/%d %.6f %.6f %.6f %.4f %.4f\n',...
         ep, epoch, loss_train(ep), loss_test(ep), P_up_est(ep), r_train, r_test);
     

end


end