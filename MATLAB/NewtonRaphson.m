function [xkpfin, storage] = NewtonRaphson(x0,func,dfunc, tol)

% INPUTS
% x0 - the initial guess
% func - the anonymous function that needs to be solved
% dfunc - the anonymous function derivative wrt the varibale be solved for
% tol - the error tolerances between xk+1 and xk for stopping


%OUTPUTS
% xkpfin is the final output
% storage stores the initial values, updated value, and iteration

err = 1;
count = 0;

while abs(err) > tol

    count = count + 1;
    xkp1 = x0 - func(x0)/dfunc(x0);
    storage(count,:) = [x0, xkp1, count];
    err = xkp1 - x0;

    x0 = xkp1;

end

xkpfin = xkp1;