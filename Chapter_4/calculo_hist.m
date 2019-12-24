inputVector11 = x;
s11{1} = double(inputVector11);
         
for r=1:size(layers,2)-2 
    e11{r} = W{r}*s11{r};
    s11{r+1} = [ones(1,60000); activation_fun(e11{r})];
         
end
outputActualInput11 = W{size(layers,2)-1}*s11{size(layers,2)-1};
for nn=1:60000
    outputVector11(:,nn) = output_fun(outputActualInput11(:,nn));   
end

cc1 = mean(abs(outputVector11 - y));
histogram(cc1,'Normalization','pdf')