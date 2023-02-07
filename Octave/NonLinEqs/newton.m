function y = Ne(x)
	func =(x^2 - 2)*e^(-x^2) %the function you want to evaluate
	funcderiv =-2*e^(-x^2)*x^3+6*EXP(x^2)*x %function derivative
	y = x - func/funcderiv
end
int i;
float newx, oldx, err
err = 10;
newx=1;
while err>0:
	newx = Ne(newx);
	err =abs(newx-oldx)
	oldx = newx
	i++;
end
