%% DK-iteration

% Importiamo il sistema dal file "Modello schema a blocchi.mat"

Delta = load("Modello schema a blocchi.mat","Delta").Delta;
G_attuata_nom = load("Modello schema a blocchi.mat", ...
                                            "G_attuata_nom").G_attuata_nom;
wi = load("Pesi incertezza.mat","W").W;
Wi = [wi,0;
      0,wi];

s = tf('s');

% Pesi sulla funzione di sensitivity
wb_p1 = 0.01;
wp1 = 0.55*wb_p1/(s+wb_p1);

wp2 = 0.1;

Wp = [wp1, 0;
      0, wp2];
Wpss = ss(Wp);
%% Plot pesi sulla funzione di sensitività

% Pesi sullo sforzo di controllo u
wu = 0.01*s/(s+0.1);
Wu = [wu, 0;
      0, wu];
Wuss = ss(Wu);

%% Plot pesi sullo sforzo di controllo u

% Pesi sulla funzione a ciclo chiuso T
Wt = 0.5*s/(s+100)*eye(2);
Wtss = ss(Wt);

%% Plot pesi sulla funzione a ciclo chiuso T


% creo la variabili per sysic
systemnames = 'G_attuata_nom Wpss Wuss Wtss Wi';
inputvar = '[u_delta(2);r(2);d(2);n(2);u(2)]'; % sarebbe [w u]'
input_to_G_attuata_nom = '[u]';
input_to_Wtss = '[G_attuata_nom + u_delta]';
input_to_Wi = '[u]';
input_to_Wpss = '[r-G_attuata_nom-d-u_delta]';
input_to_Wuss = '[u]';
outputvar = '[Wi;Wpss;Wuss;Wtss;r-G_attuata_nom-u_delta-d-n]'; % sarebbe [z v]' -- Prende così l'uscita di G + il resto

sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
P_Delta = lft(Delta,P);

opts_musyn = musynOptions('MaxIter',30,'TargetPerf',0.7);
[K_musyn,CLperf_musyn,info_musyn] = musyn(P_Delta,2,2,opts_musyn);

N_mu = lft(P,K_musyn);
BlockStructure = [2 0;6 6];
[bounds,mu_info] = mussv(ss(N_mu),BlockStructure,'o');
mu = squeeze(bounds.ResponseData);
w_mu = bounds.Frequency;

% loglog(w_mu,mu);
N_tot = lft(P_Delta, K_musyn);

%% Mu Analysis

% Stabilità
    % Nominale
[stabmarg_N_mu,wcu_N_mu,info_N_mu] = robuststab(N_mu); 
    % Robusta
[stabmarg_N_mu_p,wcu_N_mu_p,info_N_mu_p] = robuststab(N_tot);

%Prestazioni
    % Nominali
[perfmarg_N_mu,wcu_perf_N_mu,info_perf_N_mu]=robustperf(N_mu(3:end,3:end));
    % Robuste
[perfmarg_N_mu_p,wcu_perf_N_mu_p,info_perf_N_mu_p] = robustperf(N_tot);




%% Plot 
% figure('name','Valori singolari N nominale DK-iteration')
% sigma(N_mu(3:end,3:end));
% legend('N nominale');
% 
% figure('name','Valori singolari N perturbato DK-iteration')
% sigma(N_tot);
% legend('N perturbata');
% 
% figure
% plot(out.phi_plot_h_inf);
% str = 'Andamento di ${\phi}$ per modello lineare';
% title(str,'Interpreter','latex')
% xlabel('Time [s]')
% ylabel('Gradi [°]')
% 
% figure
% plot(out.phi_dot_plot_h_inf);
% str = 'Andamento di ${\dot\phi}$ per modello lineare';
% title(str,'Interpreter','latex')
% xlabel('Time [s]')
% ylabel('Gradi [°]')
%
figure
plot(out.phi_nl);
str = 'Andamento di ${\phi}$ per modello non lineare';
title(str,'Interpreter','latex')
xlabel('Time [s]')
ylabel('Gradi [°]')

figure
plot(out.phi_dot_nl);
str = 'Andamento di ${\dot\phi}$ per modello non lineare';
title(str,'Interpreter','latex')
xlabel('Time [s]')
ylabel('Gradi [°]')





