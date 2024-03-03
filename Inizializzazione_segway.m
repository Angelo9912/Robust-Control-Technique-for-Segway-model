%% Script di inizializzazione progetto CSI Brigida-Calzaretta-Massara

% definiamo i parametri del nostro sistema
g = 9.81;
syms m r a L Km Tm M
syms x y theta phi v theta_dot phi_dot x_dot y_dot theta_dot v_dot...
    theta_dotdot phi_dotdot tau_r tau_l

M_vdot = 2*(3/2*m +0.5*M);
M_theta2dot = 2*(3/2*m*a^2 + 1/6*M*a^2 + 1/2*M*L^2*sin(phi)^2);
M_phi2dot = 2*(1/2*m*r^2 + 2/3*M*L^2);
M_vdot_phi2dot = -m*r+M*L*cos(phi);
c_v = -M*L*sin(phi)*(theta_dot^2 + phi_dot^2);
c_thetadot = 2*M*L^2*sin(phi)*cos(phi)*theta_dot*phi_dot+M*L*sin(phi)...
                *v*theta_dot;
c_phidot = -M*L^2*sin(phi)*cos(phi)*theta_dot^2;
g_phidot = -M*L*g*sin(phi);

I = [M_vdot 0 M_vdot_phi2dot;
     0 M_theta2dot 0;
     M_vdot_phi2dot 0 M_phi2dot];

C = [c_v c_thetadot c_phidot]';
G = [0 0 g_phidot]';

Ing = [1/r 1/r;
       a/r -a/r;
       -1 -1];

% Equazioni delle dinamica
csi1_dot = v*cos(theta);
csi2_dot = v*sin(theta);
csi3_dot = theta_dot;
csi4_dot = phi_dot;
csi567_dot = inv(I)*(Ing*[tau_r,tau_l]' - C -G);

csi_dot = [csi1_dot; csi2_dot;csi3_dot;csi4_dot;csi567_dot];
csi = [x,y,theta,phi,v,theta_dot,phi_dot];

% Linearizzazione attorno a csi_eq
A_sym = jacobian(csi_dot,csi);
B_sym = jacobian(csi_dot,[tau_r,tau_l]);

% Modello lineare dei sensori (phi e phi_dot)
C = [0 0 0 1 0 0 0
     0 0 0 0 0 0 1];

D = zeros(2,2);

% Punto di lavoro
csi_eq = [0 0 0 0 0 0 0];

A_eq = subs(A_sym,csi,csi_eq);
B_eq = subs(B_sym,csi,csi_eq);

m = 1.22;
r = 0.13;
a = 0.496/2;
L = 0.2924;
Km = 1.08;
Tm = 0.005;
g = 9.81;
M = 54;

% Parametri incerti
u_L = ureal('u_L',L,'Percentage',[-10 +10]);
u_Km = ureal('u_Km',Km,'Percentage',[-10 +10]);
u_Tm = ureal('u_Tm',Tm,'Percentage',[-20 +20]);


A = [0 0 0 0 0 1 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 1 0;
     0 0 0 0 0 0 1;
     0 0 0 (2943*u_L*M*(m*r - u_L*M))/(100*(u_L^2*M^2+12*u_L^2*M*m+6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2)) 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 (2943*u_L*M*(M + 3*m))/(100*(u_L^2*M^2+12*u_L^2*M*m + 6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2)) 0 0 0];

B = [                                                                                                                                                         0,                                                                                                                                                              0;
                                                                                                                                                              0,                                                                                                                                                               0;
                                                                                                                                                              0,                                                                                                                                                               0;
                                                                                                                                                              0,                                                                                                                                                               0;
(4*M*u_L^2 + 3*m*r^2)/(r*(u_L^2*M^2 + 12*u_L^2*M*m + 6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2)) - (3*(m*r - u_L*M))/(u_L^2*M^2 + 12*u_L^2*M*m + 6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2), (4*M*u_L^2 + 3*m*r^2)/(r*(u_L^2*M^2 + 12*u_L^2*M*m + 6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2)) - (3*(m*r - u_L*M))/(u_L^2*M^2 + 12*u_L^2*M*m + 6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2);
                                                                                                                                    (3*a)/(r*(M*a^2 + 9*a^2*m)),                                                                                                                                    -(3*a)/(r*(M*a^2 + 9*a^2*m));
      (3*(m*r - u_L*M))/(r*(u_L^2*M^2 + 12*u_L^2*M*m + 6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2)) - (3*(M + 3*m))/(u_L^2*M^2 + 12*u_L^2*M*m + 6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2),       (3*(m*r - u_L*M))/(r*(u_L^2*M^2 + 12*u_L^2*M*m + 6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2)) - (3*(M + 3*m))/(u_L^2*M^2 + 12*u_L^2*M*m + 6*u_L*M*m*r + 3*M*m*r^2 + 6*m^2*r^2)];


% Modello "state space" del sistema perturbato (senza attuatori)
G_p = ss(A,B,C,D);

s = tf('s');

% Funzione di trasferimento dell'attuatore
Gm = u_Km/(u_Tm*s + 1);

G_tau = [Gm,0;
         0,Gm];

% Matrice di trasferimento del sistema totale
G_attuata_p = G_p*G_tau;

% Modello di incertezza additiva
E = (G_attuata_p-G_attuata_p.Nominalvalue);

E_samples = usample(E,200);

% n_points = 20;
% ord = 2;
% 
% sigma(E_samples);
% [freq,resp_db] = ginput(n_points); % pick n points
% for i = 1:n_points                    % Converts the logarithmic
%     resp(i) = 10^(resp_db(i)/20); % response to magnitude
% end % response
% sys = frd(resp,freq); % creates frd object
% W = fitmagfrd(sys,2); % fits the frequency
% % response

% Abbiamo ottenuto un peso per approssimare le incertezze
Wi = load("Pesi incertezza.mat").W;
Wi = [Wi,0;
      0, Wi];

% Modello "state space" del sistema non perturbato (senza attuatori)
SYS = ss(G_p,'min');
A_nom = SYS.A;
B_nom = SYS.B;
C_nom = SYS.C;

% Matrice di trasferimento e modello "state space" degli attuatori
% (nominale)
Gm_nom = Gm.NominalValue;
G_mss = ss(Gm_nom, 'min');


Am = G_mss.A;
Bm = G_mss.B;
Cm = G_mss.C;
Dm = G_mss.D;

A_att = Am*eye(2);
B_att = Bm*eye(2);
C_att = Cm*eye(2);


% Modello "state space" del collegamento sistema + attuatori
A_tot = [A_nom B_nom*C_att; 
         zeros(2,2) A_att];
B_tot = [zeros(2,2); B_att];

C_tot = [C_nom zeros(2,2)];
D_tot = zeros(2,2);

G_attuata_nom = ss(A_tot,B_tot,C_tot,D_tot);

% Modello delle incertezze
delta = ultidyn('delta',[1 1]);
Delta = [delta,0;0,delta];

save('Modello schema a blocchi.mat',"Delta","G_attuata_nom");
save("Matrici spazio di stato.mat","A_tot","B_tot","C_tot","D_tot");
figure
sigma(E_samples,Wi);

% 
% figure
% bodemag(G_p);
% title('Modulo del sistema perturbato senza attuatori');
% 
% figure
% bodemag(G_tau);
% title('Modulo della fdt degli attuatori');
% 
% figure
% bodemag(G_attuata_p)
% title('Modulo della fdt del sistema');
% 
% figure
% bodemag(G_attuata_nom);
% title('Modulo della fdt del sistema nominale totale');
% 
