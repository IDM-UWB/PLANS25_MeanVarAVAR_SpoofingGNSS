function [AVAR_mean,AVAR_var] = AVAR_meanVar(M,N,var_r,var_e,cov_er,Beta,T,S,R)
%
% ION 2025
% Spoofing Detection in Time Domain for Stationary and Moving Receiver
% Oliver Kost, Ondøej Straka, Ondøej Dohnal, and Jindøich Duník
%
% M interval length     % M = 1,2,...
% N number of intervals % N = fix((tau-1)/M)
%
    AVAR_mean = var_e*(2*M^2+1)/(6*M) +var_r/(M*T^2) -cov_er/(M*T) +(3*S+(3-4*Beta^M+Beta^(2*M))*R)/(M^2*T^2);
   
    AVAR_var = var_r^2*(3*N - 4)/(M^2*T^4*(N-1)^2); % var_r^2
    AVAR_var = AVAR_var + cov_er^2*(3*N - 4)/(M^2*T^2*(N-1)^2); %  cov_er^2
    AVAR_var = AVAR_var + var_e^2*(3*N + 6*N*M^2 + 9*N*M^4 - 4*M^2 - 10*M^4 - 4)/(36*M^2*(N - 1)^2); % var_e^2
    AVAR_var = AVAR_var + var_r*var_e*( 3*N + 3*N*M^2 - 2*M^2 - 4)/(3*M^2*T^2*(N - 1)^2);  % var_r*var_e
    AVAR_var = AVAR_var + -var_r*cov_er*2*(3*N - 4)/(M^2*T^3*(N - 1)^2);  % var_r*cov_er
    AVAR_var = AVAR_var + -var_e*cov_er*(3*N + 3*N*M^2 - 2*M^2 - 4)/(3*M^2*T*(N - 1)^2); % var_e*cov_er

    AVAR_var = AVAR_var + S^2*(35*N - 53 + 0^(N - 2) )/(M^4*T^4*(N - 1)^2);  % S^2
    AVAR_var = AVAR_var + var_r*S*4*( 5*N - 7)/(M^3*T^4*(N - 1)^2);  % var_r*S
    AVAR_var = AVAR_var + var_e*S*2*( 5*N + 4*N*M^2 - 2*M^2 - 7)/(3*M^3*T^2*(N - 1)^2);  % var_e*S
    AVAR_var = AVAR_var + -cov_er*S*4*( 5*N - 7)/(M^3*T^3*(N - 1)^2);  % cov_er*S

    AVAR_var = AVAR_var + -var_r*R*2*(Beta^M - 1)*(10*N + 8*Beta^M + N*Beta^(2*M) - 2*Beta^(2*M) - 5*N*Beta^M - 14)/(M^3*T^4*(N - 1)^2); %  var_r*R
    AVAR_var = AVAR_var + cov_er*R*2*(Beta^M - 1)*(10*N + 8*Beta^M + N*Beta^(2*M) - 2*Beta^(2*M) - 5*N*Beta^M - 14)/(M^3*T^3*(N - 1)^2);  % cov_er*R
    AVAR_var = AVAR_var + var_e*R*(Beta^M - 1)*(2*M^2*Beta^M - 8*Beta^M - N*Beta^(2*M) - 10*N - 8*N*M^2 + 2*Beta^(2*M) - 2*M^2*Beta^(2*M) + 4*M^2 + 5*N*Beta^M + N*M^2*Beta^M + N*M^2*Beta^(2*M) + 14)/(3*M^3*T^2*(N - 1)^2);  % var_e*R
    AVAR_var = AVAR_var + S*R*(Beta^M - 1)*(12*Beta^M - 70*N + 2*0^(N - 2)*(Beta^M - 1)^3 + Beta^M*(6*N - 18) + Beta^M*(36*N - 72) + Beta^(3*M)*(2*N - 6) - Beta^(2*M)*(6*N - 18) - Beta^(2*M)*(8*N - 16) + 106)/(M^4*T^4*(N - 1)^2); % S*R

    if N==2 
        F_t = [ 2, -12, 18];
    elseif N==3
        F_t = [ 1, -6, 21, -48, 52];
    else
        F_t=[1, -6, 17, -32, (48+16*(0:2*N-6)).*(-1).^((0:2*N-6))];
        F_t(end)=-53 + 35*N;
        F_t(end-1)= -42*N + 78;
        F_t(end-2)=-85 + 35*N;
    end
    AVAR_var = AVAR_var + R^2*(Beta^M - 1)^2*Beta.^(M*(2*N-2:-1:0))*F_t'/(T^4*M^4*(N-1)^2); % R^2  
end