function [y] = Softmax(z)

for j = 1:size(z,1)
    y(j,1) = exp(z(j))./sum(exp(z));
end

end

