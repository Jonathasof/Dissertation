function K = km_kernel(X1,X2,type,width)

switch type
	case 'gauss'	% Gaussian kernel
		sgm = 1/(2*width^2);	
		
		dim1 = size(X1,2);
 		dim2 = size(X2,2);
		
		norms1 = sum(X1.^2,1);
		norms2 = sum(X2.^2,1);
		
		mat1 = repmat(norms1',1,dim2);
		mat2 = repmat(norms2,dim1,1);
		
		distmat = mat1 + mat2 - 2*X1'*X2;
		K = exp(-distmat*sgm);
		
    
	case 'poly'	% polynomial kernel
		p = width(1);	
		c = width(2);	
		
		K = (X1'*X2 + c).^p;
		
	case 'linear' % linear kernel
		K = X1'*X2;
		

end