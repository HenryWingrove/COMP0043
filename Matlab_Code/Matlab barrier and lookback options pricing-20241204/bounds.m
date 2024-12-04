function [low,up] = bounds(parameters)

% Parameters
maxnummoments = 10;
factor_lowerbarrier = 5;
factor_upperbarrier = 5;
tol = 1e-8;

% Levy moments
[cumul1,cumul2] = cumulants(parameters);

% Lower bound
low = cumul1 - factor_lowerbarrier*sqrt(cumul2);
bound = BoundLowerTailLevy(-low,maxnummoments,parameters);
%fprintf('lower limit %.12f bound %.12f\n',low, bound);
while bound > tol
    factor_lowerbarrier = factor_lowerbarrier + 1;
    low = cumul1 - factor_lowerbarrier*sqrt(cumul2);
    bound = BoundLowerTailLevy(-low,maxnummoments,parameters);
    %fprintf('lower bound %f %.12f %.12f \n',factor_lowerbarrier,low,bound);
end
	
% Upper bound
up = cumul1 + factor_upperbarrier*sqrt(cumul2);
bound = BoundUpperTailLevy(up,maxnummoments,parameters);
%fprintf('upper limit %.12f bound %.12f\n',up,bound);
while bound > tol
    factor_upperbarrier = factor_upperbarrier+1;
    up = cumul1 + factor_upperbarrier*sqrt(cumul2);
    bound = BoundUpperTailLevy(up,maxnummoments,parameters);
    %fprintf('upper bound %f %.12f %.12f \n',factor_upperbarrier,up,bound);
end

function minup = BoundUpperTailLevy(x,maxnummoments,parameters);

minup = 1;
for j = 1:maxnummoments
    bound = real(charfunction(-j*i,parameters))*exp(-x*j);
    minup = min(minup,bound);
end

function minlow = BoundLowerTailLevy(x,maxnummoments,parameters);

minlow = 1;
for j = 1:maxnummoments
    bound = real(charfunction(j*i,parameters))*exp(-x*j);
    minlow = min(minlow,bound);
end
