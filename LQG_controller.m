%% Controllo LQG

% Importiamo le matrici della dinamica del sistema
A_tot = load("Matrici spazio di stato.mat","A_tot").A_tot;
B_tot = load("Matrici spazio di stato.mat","B_tot").B_tot;
C_tot = load("Matrici spazio di stato.mat","C_tot").C_tot;
D_tot = load("Matrici spazio di stato.mat","D_tot").D_tot;

%% Controllo senza integratore
% Definiamo i rumori di misura e i disturbi in ingresso

Qk = 0.01*diag([1 1]);          % Covarianza dei disturbi in ingresso
Rk = (0.035)^2*diag([1 1]);     % Covarianza dei rumori di misura

% Definiamo i pesi del funzionale di costo (LQR)
Qr = diag([(20*pi/180)^-2 (20*pi/180)^-2 1 1]);
Rr = (0.1)^-2*eye(2);

% Troviamo un controllore ottimo LQR
K = -lqr(A_tot,B_tot,Qr,Rr);

% Sistema in spazio di stato
G_tot = ss(A_tot,B_tot,C_tot,D_tot);

% Troviamo il guadagno del filtro di Kalman
[kalmf,L,P] = kalman(G_tot,Qk,Rk);

%% Controllo con Integratore LQI
p = 2;
n = rank(ctrb(A_tot,B_tot));
Qri = [Qr zeros(4,2);
       zeros(2,4) diag([1 1])];
Rri = Rr;
Aug = [A_tot,B_tot;
       -C_tot,zeros(2,2)];

is_possible_LQI = (rank(Aug) == n + p); % Se è true allora è possibile
                                        % fare la sintesi LQI

% non è possibile realizzare il controllo LQG con integratore in quanto
% passando al sistema in variabili aumentate (phi,phi_dot e integrali 
% del segnale di errore), otteniamo che esso non è
% stabilizzabile e quindi non è possibile sintetizzare il controllore
% (questo vale sia per il caso lineare che per quello non lineare)

figure
plot(out.phi_LQG_NL);
title('Andamento di Phi per sistema non lineare');

figure
plot(out.phi_dot_LQG_NL)
title('Andamento di Phi_dot per sistema non lineare');

% figure
% plot(out.phi_LQG_LIN);
% title('Andamento di Phi per sistema lineare');
% 
% figure
% plot(out.phi_dot_LQG_LIN)
% title('Andamento di Phi_dot per sistema lineare');