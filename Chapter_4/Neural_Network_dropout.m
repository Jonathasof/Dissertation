function [W, loss_train, loss_test] = Neural_Network_dropout(x,y,layers,learning_rate,...
                                        activation_fun,der_activation_fun,...
                                        output_fun,der_output_fun,num,loss_fun,epoch,K, x_test, y_test, d2)



c = 0;

% The number of training vectors.
trainingSetSize = size(x, 2);
d1 = 1;
%d2 = 0.8;
x = [ones(1,size(x,2));x];

% Initialize the weights for the hidden layer and the output layer.  
for i=1:size(layers,2)-1
    W{i} = randn(layers(i+1),layers(i));
    delta_w{i} = zeros(layers(i+1),layers(i)+1);
%     delta_w_antigo{i} = zeros(layers(i+1),layers(i)+1);
    W{i} = [0.01*ones(layers(i+1),1) W{i}]; %come�ando com bias zerado
    if i ==1
        W{i} = W{i}./sqrt(d1*layers(i)/(2));
    else
        W{i} = W{i}./sqrt(d2*layers(i)/(2));
    end
end 

one2 = ones(layers(end)+1,1);
iter = floor(trainingSetSize/K);

for ep = 1: epoch
%     H    = waitbar(0, 'Please wait...');
%     c = 0;
    for t = 1: iter
%         if ~ishghandle(H)
%            error('User closed waitbar.');
%         end
%         c = c + 1;     
%         H = waitbar(c / iter, H, sprintf('%d of %d', c, iter));
        drop{1} = find([1; binornd(1,d1,layers(1),1)]);%find(one1);
        drop{2} = find([1; binornd(1,d2,layers(1+1),1)]);
        drop{3} = find([1; binornd(1,d2,layers(2+1),1)]);
        drop{4} = find([1; binornd(1,d2,layers(3+1),1)]);
%         drop{5} = find([1; binornd(1,d2,layers(4+1),1)]);
        drop{size(layers,2)} = find(one2);
        
         % Select which input vector to train on.
         n = randperm(size(x,2),K);%floor(rand(1,K)*trainingSetSize + 1);
            
         % Propagate the input vector through the network.
         inputVector = x(:, n);
         s{1} = inputVector(drop{1},:);
         for r=1:size(layers,2)-2 
             z{r} = W{r}(drop{r+1}(2:end)-1,drop{r})*s{r};
             s{r+1} = [ones(1,K); activation_fun(z{r})];
         end
         
         outputActualInput = W{size(layers,2)-1}(drop{end}(2:end)-1,drop{end-1})*s{size(layers,2)-1};
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
         for r=size(layers,2)-2:-1:1
             er{r} = W{r+1}(drop{r+2}(2:end)-1,drop{r+1}(2:end))'*delta{r+1};
             delta{r} = der_activation_fun(z{r}).*(er{r});  
         end 
        

         for i = 1:size(layers,2)-2
              delta_w{i+1}(drop{i+2}(2:end)-1,drop{i+1}) =  - (learning_rate/K)*delta{i+1}*s{i+1}';
              W{i+1}(drop{i+2}(2:end)-1,drop{i+1}) = W{i+1}(drop{i+2}(2:end)-1,drop{i+1}) + delta_w{i+1}(drop{i+2}(2:end)-1,drop{i+1});
         end
         delta_w{1}(drop{2}(2:end)-1,drop{1}) = - (learning_rate/K).*delta{1}*s{1}';
         W{1}(drop{2}(2:end)-1,drop{1}) = W{1}(drop{2}(2:end)-1,drop{1}) + delta_w{1}(drop{2}(2:end)-1,drop{1});


     end
        % Calculate the loss and accuracy for plotting.
         s{1} = x;
         s_test{1} = x_test;
         s_test{1} = [ones(1,size(s_test{1},2));s_test{1}];

         targetVector = y'; 
         targetVector_test = y_test';
         for r=1:size(layers,2)-2 
             if r == 1
                 z{r} = d1*W{r}*s{r};
                 z_test{r} = d1*W{r}*s_test{r};
             else
                 z{r} = d2*W{r}*s{r};
                 z_test{r} = d2*W{r}*s_test{r};
             end
             s{r+1} = [ones(1,size(x,2)); activation_fun(z{r})];
             s_test{r+1} = [ones(1,size(x_test,2)); activation_fun(z_test{r})];
         end
         outputActualInput = d2*(W{size(layers,2)-1})*s{size(layers,2)-1};
         for n=1:size(x,2)
              outputVector_train(:,n) = output_fun(outputActualInput(:,n));   
         end
         
         outputActualInput_test = (d2*W{size(layers,2)-1})*s_test{size(layers,2)-1};
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
% close(H)
fprintf('Epoch %d/%d; loss train: %.6f; loss test: %.6f, %.4f; %.4f \n',...
         ep, epoch, loss_train(ep), loss_test(ep), r_train, r_test);
    
end


end



