function [y] = der_Softmax(z)

for j = 1:size(z,1)
    y(j,1) = [exp(z(j))*(sum(exp(z))-exp(z(j)))]./(sum(exp(z)))^2;
end

end

