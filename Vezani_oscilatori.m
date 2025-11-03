%% Vezani jednostavni hamonijski oscilatori
% Usporedba eksplicitne Eulerove metode i simplektičke Eulerove metode

clear; clc; close all;

% --- Parametri ---
m = 1;           % masa
k = 1;           % konstanta opruge
kc = 0.2;        % konstanta zajedničke opruge
dt = 0.02;       % vremenski korak
T = 200;         % ukupno vrijeme
t = 0:dt:T;      
N = length(t);

% --- Početni uvjeti ---
x1 = zeros(1,N); x2 = zeros(1,N);
v1 = zeros(1,N); v2 = zeros(1,N);
x1(1) = 1;       % pomak prve mase
x2(1) = 0;       % druga miruje
v1(1) = 0;       v2(1) = 0;

x1s = x1; x2s = x2;  
v1s = v1; v2s = v2;

% --- Funkije sile ---
F1 = @(x1,x2) -k*x1 - kc*(x1 - x2);
F2 = @(x1,x2) -k*x2 - kc*(x2 - x1);

% --- Integracija eksplicitna ---
for i = 1:N-1
    a1 = F1(x1(i),x2(i))/m;
    a2 = F2(x1(i),x2(i))/m;

    x1(i+1) = x1(i) + dt*v1(i);
    x2(i+1) = x2(i) + dt*v2(i);
    v1(i+1) = v1(i) + dt*a1;
    v2(i+1) = v2(i) + dt*a2;
end

% --- Integracija simplektička ---
for i = 1:N-1
    a1 = F1(x1s(i),x2s(i))/m;
    a2 = F2(x1s(i),x2s(i))/m;
    
    % Prvo brzine
    v1s(i+1) = v1s(i) + dt*a1;
    v2s(i+1) = v2s(i) + dt*a2;
    
    % Onda položaj s osvježenim brzinama
    x1s(i+1) = x1s(i) + dt*v1s(i+1);
    x2s(i+1) = x2s(i) + dt*v2s(i+1);
end

% --- Energija ---
E_euler = 0.5*m*(v1.^2 + v2.^2) + 0.5*k*(x1.^2 + x2.^2) + 0.5*kc*(x1 - x2).^2;
E_symp = 0.5*m*(v1s.^2 + v2s.^2) + 0.5*k*(x1s.^2 + x2s.^2) + 0.5*kc*(x1s - x2s).^2;

% --- Plot ---
figure('Name','Coupled SHOs','Position',[100 100 1200 700])

subplot(2,2,1)
plot(t,x1,'r',t,x2,'b')
xlabel('Vrijeme'); ylabel('Pomak')

subplot(2,2,2)
plot(t,x1s,'r',t,x2s,'b')
xlabel('Vrijeme'); ylabel('Pomak')

subplot(2,2,3)
plot(t,E_euler,'LineWidth',1.5)
xlabel('Vrijeme'); ylabel('Energija')

subplot(2,2,4)
plot(t,E_symp,'LineWidth',1.5)
xlabel('Vrijeme'); ylabel('Energija')




