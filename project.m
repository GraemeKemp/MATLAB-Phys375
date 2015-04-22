format long

%lists the constants used in the equtions solved in (the other file??)
G = 6.67384e-11; %gravitational constant
h_bar = 1.054571726e-34; %reduced planck constant
m_e = 9.10938291e-31; %mass of an electron
m_p = 1.67262178e-27; %mass of a proton
stefan_boltz = 5.670e-8; %sigma
c = 2.99792458e8; %speed of light
a = (4*stefan_boltz) / c; %a thing
k = 1.3806488e-23; %boltzmann constant
Gamma = (5/3);%assuming some sort of polytropic index

X = 0.7; %1 - (2*10^-5);%distribution of masses by element
X_CNO = 0.03*X;
Y = 0.28;
Z = 1-X-Y; 
mu = (2*X + 0.75*Y + 0.5*Z)^(-1);

%T_c = 35e7;
M_sun = 1.989e30;%test values? Known Solar constants
L_sun = 3.846e26;
R_sun = 6.95800e8;

eps_abs = 1e-5;
eps_step = 1e-5;
%new part:27-29, runs 400 stars - determines the central temperatures we
%are solving the star for in order to create a main sequence

% A = 0.2e6:5e5:10e6;
% B = 10e6:2.5e5:100e6;%reset to 400 for full run(2.5e5 for B)
% T_c = cat(2,A,B);

% T_c = B;%FOR SHORT RUNS
T_c = 27e6; %was 35e7 - think its high  -try new T_c for different mass

length = size(T_c, 2);
L_star = ones(1, length);
M_star = ones(1, length);
R_star = ones(1, length);
T_star = ones(1, length);
%L_star = 2.6437e3*L_sun;
%R_star = 3.1341*R_sun;
%M_star = 6.9548*M_sun;
%T_star = 7.1020e3;
%% Mass = ones(1,length);

for j=1:length 
    rho_c_min = 300;
    rho_c_max = 500000;
    function_rho_c_new = 100000;
    i=0;
    while ((abs(real(function_rho_c_new)) > eps_abs) && (i<200))
        rho_c_new = real((rho_c_min + rho_c_max)/2);
        [function_rho_c_min, ~, ~, ~, ~, ~, ~, ~, ~, ~] = getErrorInDensity(rho_c_min,T_c(1, j));
        [function_rho_c_max, ~, ~, ~, ~, ~, ~, ~, ~, ~] = getErrorInDensity(rho_c_max,T_c(1, j));
        [function_rho_c_new, R_star_new, T_star_new, L_star_new, M_star_new, R, Rho, Temp, Mass, Lum] = getErrorInDensity(rho_c_new,T_c(1, j));
%         if(real(function_rho_c_min)*real(function_rho_c_max) > 0)
%            rho_c_new = 0;
%            R_star_new = 0;
%            T_star_new = 0;
%            L_star_new = 0;
%            M_star_new = 0;
%            R = 0;
%            Rho = 0;
%            Temp = 0;
%            Mass = 0;
%            Lum = 0;
%            break;
%         end
        if (real(function_rho_c_new) == 0)%error function loop; f(rho_c)
           break;
        elseif ( function_rho_c_new > 0 )%bisection function??
           rho_c_max = rho_c_new;
        else
           rho_c_min = rho_c_new;
        end
        if ( rho_c_max - rho_c_min < eps_step )
            if ( abs( real(function_rho_c_min) ) < abs( real(function_rho_c_max) ) && abs( real(function_rho_c_min) ) < eps_abs )
                rho_c_new = rho_c_min;
                break;
            elseif ( abs( real(function_rho_c_max) ) < eps_abs )
                rho_c_new = rho_c_max;
                break;
            end
        end
        i=i+1;
    end
    i = 0;
%LUMINOSITY
    if((L_star_new == 0) || (abs(real(function_rho_c_new)) > 1000))%new, double line is OR, why 7.5???
        L_star(1, j) = nan;
    else   
        L_star(1, j) = real(L_star_new/L_sun);
    end;
%TEMPERATURE
    if((T_star_new == 0) || (abs(real(function_rho_c_new)) > 1000))%new
        T_star(1, j) = nan;
    else   
        T_star(1, j) = real(T_star_new);
    end;
%RADIUS
    if((R_star_new == 0) || (abs(real(function_rho_c_new)) > 1000))%new
        R_star(1, j) = nan;
    else   
        R_star(1, j) = real(R_star_new/R_sun);
    end; 
% %MASS    
    if((M_star_new == 0) || (abs(real(function_rho_c_new)) > 1000))%new
        M_star(1, j) = nan;
    else   
        M_star(1, j) = real(M_star_new/M_sun);
    end;
% MASS (R)
%     if((Mass == 0) || (abs(real(function_rho_c_new)) > 7.5))%new
%         Mass(1, j) = nan;
%     else   
%         Mass(1, j) = real(Mass);
%     end; 
% RHO
%     if((Rho == 0) || (abs(real(function_rho_c_new)) > 7.5))%new
%         Rho(1, j) = nan;
%     else   
%         Rho(1, j) = real(M_star_new/M_sun);
%     end; 
end

% plot((R./R_star_new),Rho./rho_c_new);
% title('Rho')
% figure();
% plot((R./R_star_new), Temp./T_c);
% title('Temp');
% figure();
% plot((R./R_star_new), Mass./M_star_new);
% title('Mass');
% figure();
% plot((R./R_star_new), Lum./L_star_new);
% title('Lum');
% figure();

%PRESSURE
% [cols_size row_size] = size(R);
% Pressure = ones(cols_size,1).*((((3*pi^2)^(2/3)) / 5)*((h_bar^2) / m_e).* (Rho./m_p).^(5/3) + Rho.*(k.*Temp./ mu*m_p) + (1/3)*a.*Temp.^4);
% figure();
% plot(R./R_star_new, Pressure);
% title('Pres');

%dL/dR
% [cols_size row_size] = size(R);
% Epp = 1.07*10^(-7).*(Rho./ 10^5)*X^2.*(Temp./ 10^6).^4;
% Ecno = 8.24*10^(-26).*(Rho./ 10^5)*X*X_CNO.*(Temp./ 10^6).^(19.9);
% E = Epp + Ecno;
% dLdR = ones(cols_size,1).*(4.*pi.*R.^2.*(Rho).*E);
% figure();
% plot((R./R_star_new),(dLdR./(L_star_new./R_star_new)));
% title('dLdR');

%KAPPA
% [cols_size row_size] = size(R);
% Kappa_es = ones(cols_size,1).*(0.02*(1+X));
% Kappa_ff = (1.0.*10^24).*(Z+0.0001).*(Rho.^(0.7)).*(Temp.^(-3.5));
% Kappa_H = (2.5.*10.^(-32)).*(Z/0.02).*(Rho.^(0.5)).*(Temp.^(9));
% Kappa = ((1./Kappa_H) + (1./max(Kappa_es, Kappa_ff))).^(-1);
% hold all;
% plot(R./R_star_new, log10(Kappa_es));
% plot(R./R_star_new, log10(Kappa_ff));
% plot(R./R_star_new, log10(Kappa_H));
% plot(R./R_star_new, log10(Kappa));
% hold off;
% title('Kappa');

%CONVECTION ZONES
% [cols_size row_size] = size(R);
% Kappa_es = ones(cols_size,1).*(0.02*(1+X));
% Kappa_ff = (1.0.*10^24).*(Z+0.0001).*(Rho.^(0.7)).*(Temp.^(-3.5));
% Kappa_H = (2.5.*10.^(-32)).*(Z/0.02).*(Rho.^(0.5)).*(Temp.^(9));
% Kappa = ((1./Kappa_H) + (1./max(Kappa_es, Kappa_ff))).^(-1);
% Pressure =((((3*pi^2)^(2/3)) / 5)*((h_bar^2) / m_e).* (Rho./m_p).^(5/3) + Rho.*(k.*Temp./ mu*m_p) + (1/3)*a.*Temp.^4);
% dT_L = ones(cols_size,1).*(3*Kappa.*Rho.*Lum) ./ ((16*pi*a*c.*Temp.^3).*(R.^2));
% dT_R = (1 - (1/Gamma)).*(Temp./ Pressure).*((G*Mass.*Rho) ./ R.^2);
% dTdR = -1*min(dT_L, dT_R);
% figure();
% hold on;
% plot(R./R_star_new,dT_L./dTdR);%radiadtive
% plot((R./R_star_new),dT_R./dTdR);%convective
% % plot(R./R_star_new,dTdR);%temperature gradient
% hold off;
% title ('Low Mass - Convection Zones');

%ENERGY GENERATION
%[cols_size row_size] = size(R);
%Epp = ones(cols_size,1).*1.07*10^(-7).*(Rho./ 10^5)*X^2.*(Temp./ 10^6).^4;
%Ecno = 8.24*10^(-26).*(Rho./ 10^5)*X*X_CNO.*(Temp./ 10^6).^(19.9);
%E = Epp + Ecno;
%hold on;
%plot(R./R_star_new, Epp);
%plot(R./R_star_new, Ecno);
%plot(R./R_star_new, E);
%hold off;
%title('Energy Generation - High Mass Star');


plot(log10(T_star), log10(Luminosity), '*b');
set(gca,'xdir','reverse');
xlim([3 4.15])
ylim([-6 6])
xlim([3 4.1])
ylim([-6 5])
title('Main Sequence of Stars')
xlabel('Log Base 10 of Temperature')
ylabel('Log Base 10 of Luminosity Divided by the Luminosity of the Sun')
