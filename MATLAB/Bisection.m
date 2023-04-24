function [xm, count] = Bisection(range,func, tol)

% INPUTS:
% range - range of values that the zolution is allowed to exist in (2x1)
% func - annonymous func handle to be evaluated
% tol - tolerance of the final solution (1x1)

% OUTPUTS
% Xkpfin - the final updated position
% storage - keeps track of the possible solutions


test = 0;
x = range(1):.001:range(2);
f = func(x);
negind = f<=test;      % Creates logical vec checking for negative f(x)

if sum(negind) == 0
    test =tol;
    negind = f<=test;
end

posind = f>0;       % Creates logical vec checking for positive f(x)
xfneg = x(negind);  % Returns all x that give negative f(x);
xfpos = x(posind);  % Returns all x that give positve f(x);
xn = xfneg(1);
xp = xfpos(1);
fm = 1;
count = 0;

while abs(fm) > tol
    count = count+1;
    fn = func(xn);
    fp = func(xp);
    xm = .5*(xn+xp);
    fm = func(xm);

    if fm*fn <= test;
        xp =xm;
    else
        xn = xm;
    end
end
