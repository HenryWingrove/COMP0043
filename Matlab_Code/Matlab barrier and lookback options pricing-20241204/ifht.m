% Fast Hilbert transform: Hilbert transform trough fast Fourier transform
% Guido Germano, 2013
function iHF = ifht(F,type,P)

% Setup
[M N] = size(F); % Dimension parameters: number of equations and of grid points
Q = N + P; % Number of grid points after zero padding

% Define the auxiliary vector
persistent vec type_old; % Declare static variables
if sum(size(vec) ~= [M Q]) || ~strcmp(type_old,type) % Initialise static variables
    type_old = type;
    if     strcmp(type,'sgn1') % Right-continuous sign function
        vec = [-ones(M,Q/2) ones(M,Q/2)];
    elseif strcmp(type,'sgn0') % Symmetric sign function
        vec = [zeros(M,1) -ones(M,Q/2-1) zeros(M,1) ones(M,Q/2-1)];
    elseif strcmp(type,'sgn2') % Left-continuous sign function
        vec = [ones(M,1) -ones(M,Q/2) ones(M,Q/2-1)];
    elseif strcmp(type,'sinc') % Sinc method
        % 1. L. Feng, V. Linetsky, Math. Finance 18, 337-384 (2008):
        %    Eq. (6.39), p. 365.
        % 2. F. Stenger, Handbook of Sinc Numerical Methods, Chapman&Hall/CRC,
        %    Boca Raton, 2011: Eq. (1.5.86), p. 107.
        % 3. F. Stenger, Numerical Methods Based on Sinc and Analytic Functions,
        %    Springer, New York, 1993.
        % 4. J. A. C. Weideman, Comp. Math. 64 (210), 745-762 (1995):
        %    Kress and Martensen method.
        t = (1-(-1).^(-Q/2:Q/2-1))./(pi*(-Q/2:Q/2-1));
        t(Q/2+1) = 0;
        vec = repmat(imag(fft(ifftshift(t))),M,1);
    end
end

% Compute the Hilbert transform times the imaginary unit
f = ifft(F,Q,2); % Optional padding with P trailing zeros to length Q = N + P
iHF = fft(vec.*f,[],2);
iHF = iHF(:,1:N); % iHF = -iHF(:,1:N) if fft and ifft are swapped in the previous two lines

% for i = 1:M
%     figure, hold on
%     plot(real(F(i,1:N)),'r')
%     plot(imag(F(1,1:N)),':r')
%     plot(real(-1i*iHF(i,:)),'b')
%     plot(imag(-1i*iHF(i,:)),':b')
%     legend('Re f','Im f','Re Hf','Im Hf')
%     title('FHT')
% end