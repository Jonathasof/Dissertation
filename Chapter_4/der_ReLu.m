function [y] = der_ReLu(x)

y = min(max(0,ceil(x)),1);

end

