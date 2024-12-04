function [alpha] = set_alpha(param,call,d,u,b)

%-- Feng-Linetsky, page 351
% P_type = 0 vanilla call payoff, 1 truncated payoff
if -b < d && u == b % down-and-out option
    if call == 1 %DOC
        P_type = 0;
    else %DOP
        P_type = 1;
    end
elseif -b == d && u < b % up-and-out option
    if call == 1 %UOC
        P_type = 1;
    else %UOP
        P_type = 0;
    end
else % double-barrier option
    P_type = 1;
end

%-- Feng-Linetsky, page 364
lm=param.lambdam;
lp=param.lambdap;
if lm==0 || lp==0 %\pm\infty
    if P_type == 0
        % vanilla call payoff case
        alpha=-10; 
    else 
        % truncated payoff case
        alpha=0; 
    end
else
    if P_type == 0
        % vanilla call payoff case
        alpha=(lm-1)/2; 
    else 
        % truncated payoff case
        alpha=(lp+lm)/2; 
    end
end