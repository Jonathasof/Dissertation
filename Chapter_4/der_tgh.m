function [y]=der_tgh(x)
y = 1/2*(1-tgh(x/2).^2);
end


