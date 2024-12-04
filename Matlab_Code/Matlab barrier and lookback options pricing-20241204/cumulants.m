function [cumul1,cumul2] = cumulants(parameters)

dt = parameters.dt;
rf = parameters.rf;

switch parameters.distr

    case 1 % Normal

        s = parameters.s;

        cumul1 = dt*(rf-0.5*s^2);
        cumul2 = dt*s^2;

    case 2 % NIG, Schoutens page 59

        alpha = parameters.alpha;
        beta = parameters.beta;
        delta = parameters.delta/dt;

        cumul1 = dt*(rf+delta*((beta-alpha^2+beta^2)/sqrt(alpha^2-beta^2) ...
                 + sqrt(alpha^2-(beta+1)^2)));
        cumul2 = dt*delta*alpha^2/(alpha^2-beta^2)^1.5;

    case 3 % VG, Schoutens page 57

        nu = parameters.nu*dt;
        theta = parameters.theta/dt;
        s = parameters.s/sqrt(dt);

        cumul1 = dt*(rf+theta-log((1-nu*(0.5*s^2+theta))^(-1/nu)));
        cumul2 = dt*(s^2+nu*theta^2);

     case 4 % Meixner, Schoutens page 62

        alpha = parameters.alpha;
        beta = parameters.beta;
        delta = parameters.delta/dt;

        cumul1 = dt*(rf-log((cos(beta/2)*sec((alpha+beta)/2))^2*delta) ...
                 + alpha*delta*tan(beta/2));
        cumul2 = dt*alpha^2*delta*(sec(beta/2)^2)/2;

     case 5 % CGMY, Schoutens page 60

        C = parameters.C/dt;
        G = parameters.G;
        M = parameters.M;
        Y = parameters.Y;

        cumul1 = dt*(G*M*rf - C*(-G^(Y+1)*M - G^Y*M*Y + G*(((G+1)^Y ...
                 + (M-1)^Y)*M - M^(Y+1) + M^Y*Y))*gamma(-Y))/(G*M);
        cumul2 = dt*C*(G^Y*M^2+G^2*M^Y)*(Y-1)*Y*gamma(-Y)/(G*M)^2;

    case 6 % Kou double exponential, Fusai and Roncoroni page 53

        s = parameters.s/sqrt(dt);
        lambda = parameters.lambda/dt;
        pigr = parameters.pigr;
        eta1 = parameters.eta1;
        eta2 = parameters.eta2;

        cumul1 = dt*(eta2*lambda*pigr + eta1*(lambda*(pigr-1)+eta2*rf) ...
                 - eta1*eta2*((lambda*(1+eta1*(pigr-1) + eta2*pigr)) ...
                 /((eta1-1)*(eta2+1)) + 0.5*s^2))/(eta1*eta2);
        cumul2 = dt*(-lambda*(2*(pigr-1)/eta2^2-2*pigr/eta1^2) + s^2);

    case 7 % Merton jump-diffusion, Fusai and Roncoroni page 53

        s = parameters.s/sqrt(dt);
        alpha = parameters.alpha;
        lambda = parameters.lambda/dt;
        delta = parameters.delta;

        cumul1 = dt*((alpha+1-exp(alpha+0.5*delta^2))*lambda + rf - 0.5*s^2);
        cumul2 = dt*(lambda*(alpha^2 + delta^2) + s^2);

    case 8 % Stable

        alpha = parameters.alpha;
        beta = parameters.beta;
        gamm = parameters.gamm;
        m = parameters.m;
        c = parameters.c;

        cumul1 = dt*(rf-gamm^alpha); % true only for alpha = 2
        cumul2 = dt*2*gamm^alpha; % true only for alpha = 2

end
