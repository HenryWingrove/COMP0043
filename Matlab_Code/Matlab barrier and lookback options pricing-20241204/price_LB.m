function price = price_LB(sol,x,lx,lk,parameters)

% Correzione: se la soluzione numerica supera 1 o non è più crescente
% (causa errori di approssimazione numerica)....
temp = max(sol); temp = min(1,temp); temp = find(sol>=temp,1); % primo punto dove avviene l'errore
% modifica la soluzione nell'intervallo [x(temp), x(end)]
% 1- sostituendola con una retta che va da A a B
a = x(temp); A = min(sol(temp),1);
b = x(end); B = 1;
sol(temp:end) = (A*(x(temp:end)-b)-B*(x(temp:end)-a))./(a-b);
% 2- oppure mandando tutto a 1
%sol(temp:end) = 1;

% Calcola il prezzo dell'opzione con un'ulteriore integrazione su N nodi
ss=linspace(lx-lk,x(end),length(x))';
ww = weights(4,length(x))*(ss(2)-ss(1));
%[ss,ww] = gauleg(lx-lk,x(end),length(x));
p = 1-interp1(x,sol,ss,'spline');
price = exp(-parameters.rf*parameters.T)*sum(ww.*(exp(lx-ss).*p));

function [x,w] = gauleg(a,b,n)
% Computes the Gauss-Legendre nodes x and weights w
% of order n for the interval (a,b)

m=fix((n+1)/2);
xm=0.5*(b+a);
xl=0.5*(b-a);
xx1=zeros(1,m); ww1=zeros(1,m);
xx2=zeros(1,m); ww2=zeros(1,m);
for i=1:m
    z=cos(pi*(i-0.25)/(n+0.5));
    z1=z+2*eps; %initialization to enter the while loop
    while (abs(z-z1)>=eps)
        p1=1.0;
        p2=0.0;
        for j=1:n
            p3=p2;
            p2=p1;
            p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        end
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;
    end
    xx1(i)=xm-xl*z;
    ww1(i)=2.0*xl/((1.0-z*z)*pp*pp);
    xx2(i)=xm+xl*z;
    ww2(i)=ww1(i);
end

if m==(n+1)/2
    x=[xx1,xx2(end-1:-1:1)];
    w=[ww1,ww2(end-1:-1:1)];
else
    x=[xx1,xx2(end:-1:1)];
    w=[ww1,ww2(end:-1:1)];
end
