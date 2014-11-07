% Assignmenet Question 3
% Nicholas Speal
% Assignment 4


clear
clc
close all

%% General System Parameters
a21 = 2;
r = 2;
A = [0 1; a21 1];
B = [0;1];
C = [1 0];
D = 0;


x0 = [1;1];
x0_est = [0;0];
t = [0:.01:15]';
r = 2*ones(size(t));

n = size(A);

%% A
disp('PART A')
des_eigs = [-1 -2];
K = place(A,B,des_eigs)
G = -inv(C*inv(A-B*K)*B)

A_cl = A-B*K;
B_cl = B*G;
sys_cl = ss(A_cl,B_cl,C,D);

[y t x] = lsim(sys_cl,r,t,x0);

DC_gain_for_part_a = y(end)/r(end)

%% Part B
disp('PART B')
des_eig_obs = [-4 -5];
L = place(A',C',des_eig_obs)'

% observer dynamics equation:

A_obs = A-L*C;
B_obs = [B L];
U_obs = [r zeros(size(t))];

observer_system = ss(A_obs,B_obs,C,D);
[y_est t x_est] = lsim(observer_system, U_obs, t, x0_est);

figure
plot(t, y_est); xlabel('time'); ylabel('y_{est}'); title('Observer Dynamics');

%% PART C
disp('PART C')
% Full state feedback system with observer
A_obs_fb = [A -B*K; L*C A-B*K-L*C];
B_obs_fb = [G*B; G*B];
C_obs_fb = [C zeros(1,2)];
D_obs_fb = 0;
sys_obs_fb = ss(A_obs_fb, B_obs_fb, C_obs_fb, D_obs_fb)

[y, t, x] = lsim(sys_obs_fb, r, t, [x0; x0_est]);

DC_gain_for_part_C = y(end)/r(end)

%control input for each variable:
U = -K*x(:,3:4)' + G*r';      %transposes introduced to make dimensions agree
%ctrlInput = B_obs_fb*U; %not sure...
ctrlInput = [1;1]*U;


figure
subplot(211)
    plot(t,x(:,1),t,x(:,3),t,ctrlInput(1,:))
    ylabel('x_1','fontsize',20),legend('x_1','x_1 - estimated','Control Input')    
subplot(212)
    plot(t,x(:,2),t,x(:,4),t,ctrlInput(2,:))
    ylabel('x_2','fontsize',20),legend('x_2','x_2 - estimated','Control Input')
    
%% PART D
disp('PART D')
a21 = 2.1;
A = [0 1; a21 1]; 

% Full state feedback system with observer
A_obs_fb = [A -B*K; L*C A-B*K-L*C];
B_obs_fb = [G*B; G*B];
C_obs_fb = [C zeros(1,2)];
D_obs_fb = 0;
sys_obs_fb = ss(A_obs_fb, B_obs_fb, C_obs_fb, D_obs_fb)

[y, t, x] = lsim(sys_obs_fb, r, t, [x0; x0_est]);

DC_gain_for_part_D = y(end)/r(end)

%control input for each variable:
U = -K*x(:,3:4)'+G*r';      %transposes introduced to make dimensions agree
ctrlInput = [1;1]*U;
%ctrlInput = B_obs_fb*U; %not sure...
    
figure
subplot(211)
    plot(t,x(:,1),t,x(:,3),t,ctrlInput(1,:))
    ylabel('x_1','fontsize',20),legend('x_1','x_1 - estimated','Control Input')    
subplot(212)
    plot(t,x(:,2),t,x(:,4),t,ctrlInput(2,:))
    ylabel('x_2','fontsize',20),legend('x_2','x_2 - estimated','Control Input')    
      
%% PART E
disp('PART E')
%restore nominal system
a21 = 2; 
A = [0 1; a21 1];

%augment with integral action as a servomechanism
A_aug = [A zeros(n,1); -C 0];
B_aug = [B; 0];
DesEigs = [-1 -2 -3];
K_total = place(A_aug,B_aug, DesEigs)
    k = K_total(1:end-1);
    ki = -K_total(end);

A_servo = A_aug-B_aug*K_total;
B_servo = [zeros(n,1);1];
C_servo = [C 0];
sys_servo = ss(A_servo, B_servo, C_servo, D);

[num,den]=ss2tf(A_servo,B_servo,C_servo,D);
closed_loop_poles_part_E = roots(den)
disp('negative real poles mean asymptotic stability')
tf(num,den)
DC_gain_for_part_E = num(end)/den(end)

%% PART F
disp('PART F')
%purturb system
a21 = 2.1;
A = [0 1; a21 1];

A_aug = [A zeros(n,1); -C 0];
B_aug = [B; 0];

%augment with integral action as a servomechanism
A_servo = A_aug-B_aug*K_total;  %use same K_total as in part E
B_servo = [zeros(n,1);1];
C_servo = [C 0];
sys_servo = ss(A_servo, B_servo, C_servo, D);

[num,den]=ss2tf(A_servo,B_servo,C_servo,D);
closed_loop_poles_part_F = roots(den)
disp('negative real poles mean asymptotic stability')
tf(num,den)
DC_gain_for_part_F = num(end)/den(end)


%% PART G
disp('PART G')


%--- PERTURBED 

% create big system with servomechanism, observer, and perturbed plant
% use K and L from parts E and B respectively

A_BigSys_pert = [A B*ki -B*k; -C 0 zeros(1,n); L*C B*ki A-B*k-L*C];
B_BigSys = [zeros(n,1); 1; zeros(n,1)];
C_BigSys = [C 0 zeros(1,n)];

sys_pert = ss(A_BigSys_pert, B_BigSys, C_BigSys, 0);



%--- NOMINAL

%recreate big system with servo, observer, and nominal plant
%use same K and L
a21 = 2.0;
A = [0 1; a21 1];

A_BigSys_Nominal = [A B*ki -B*k; -C 0 zeros(1,n); L*C B*ki A-B*k-L*C];
%B,C,D unchanged
sys_nominal = ss(A_BigSys_Nominal, B_BigSys, C_BigSys, 0);


%--- COMPARISON

%DC Gain Analysis
[num,den]=ss2tf(A_BigSys_pert, B_BigSys, C_BigSys, 0)
DC_gain_perturbed = num(end)/den(end)

[num,den]=ss2tf(A_BigSys_Nominal, B_BigSys, C_BigSys, 0)
DC_gain_nominal = num(end)/den(end)

%Simulate
[y_pert, t, x_pert] = lsim(sys_pert, r, t, [x0; 0; x0_est]);
[y_nom, t, x_nom] = lsim(sys_nominal, r, t, [x0; 0; x0_est]);


%control input for each variable:
U_nom = -K*x_nom(:,4:5)' +G*r';      %transposes introduced to make dimensions agree
U_pert = -K*x_pert(:,4:5)' +G*r';      %transposes introduced to make dimensions agree

    
figure
subplot(311)
    plot(t,x_nom(:,1),t,x_nom(:,4),t,U_nom,t,x_pert(:,1),t,x_pert(:,4),t,U_pert)
    ylabel('x_1','fontsize',20),legend('Nominal Plant: x_1','Nominal Plant: x_1 - estimated','Nominal Plant: Control Input', 'Perturbed Plant: x_1','Perturbed Plant: x_1 - estimated','Perturbed Plant: Control Input')    
subplot(312)
    plot(t,x_nom(:,2),t,x_nom(:,5),t,U_nom,t,x_pert(:,2),t,x_pert(:,5),t,U_pert)
    ylabel('x_2','fontsize',20),legend('Nominal Plant: x_2','Nominal Plant: x_2 - estimated','Nominal Plant: Control Input', 'Perturbed Plant: x_2','Perturbed Plant: x_2 - estimated','Perturbed Plant: Control Input')
subplot(313)
    plot(t,x_nom(:,3),t,x_pert(:,3))
    ylabel('Integral Error'),legend('Nominal Plant', 'Perturbed Plant')
    
    
    
    
    
    
    
    
    