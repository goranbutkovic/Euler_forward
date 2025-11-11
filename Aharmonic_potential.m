%% 2D particle in nonlinear potential: Explicit vs Symplectic Euler with energy plot

clear all; close all; clc;

% Parameters
m = 1;          % mass (kg)
kx = 1;         % spring constant x
ky = 1;         % spring constant y
dt = 0.01;      % time step
tf = 100;        % final time
t = 0:dt:tf;

% Nonlinear potential
V = @(x,y) 0.5*kx*x.^2 + 0.5*ky*y.^2 + 0.1*(x.^4 + y.^4);

% Initial conditions
x0 = 1; y0 = 0.5; 
vx0 = 0; vy0 = 0;

%% Explicit Euler
X_exp = zeros(length(t),4); % [x,y,vx,vy]
X_exp(1,:) = [x0,y0,vx0,vy0];
E_exp = zeros(length(t),1);

for n = 1:length(t)-1
    x = X_exp(n,1); y = X_exp(n,2);
    vx = X_exp(n,3); vy = X_exp(n,4);
    
    % Forces
    Fx = -kx*x - 0.4*x^3;
    Fy = -ky*y - 0.4*y^3;
    
    % Explicit Euler
    vx_new = vx + dt*Fx/m;
    vy_new = vy + dt*Fy/m;
    x_new = x + dt*vx;
    y_new = y + dt*vy;
    
    X_exp(n+1,:) = [x_new, y_new, vx_new, vy_new];
    
    % Energy
    KE = 0.5*m*(vx^2 + vy^2);
    PE = V(x,y);
    E_exp(n) = KE + PE;
end
E_exp(end) = 0.5*m*(X_exp(end,3)^2 + X_exp(end,4)^2) + V(X_exp(end,1),X_exp(end,2));

%% Symplectic Euler
X_sym = zeros(length(t),4); % [x,y,vx,vy]
X_sym(1,:) = [x0,y0,vx0,vy0];
E_sym = zeros(length(t),1);

for n = 1:length(t)-1
    x = X_sym(n,1); y = X_sym(n,2);
    vx = X_sym(n,3); vy = X_sym(n,4);
    
    % Forces
    Fx = -kx*x - 0.4*x^3;
    Fy = -ky*y - 0.4*y^3;
    
    % Symplectic Euler
    vx_new = vx + dt*Fx/m;
    vy_new = vy + dt*Fy/m;
    x_new = x + dt*vx_new;
    y_new = y + dt*vy_new;
    
    X_sym(n+1,:) = [x_new, y_new, vx_new, vy_new];
    
    % Energy
    KE = 0.5*m*(vx_new^2 + vy_new^2);
    PE = V(x_new, y_new);
    E_sym(n) = KE + PE;
end
E_sym(end) = 0.5*m*(X_sym(end,3)^2 + X_sym(end,4)^2) + V(X_sym(end,1),X_sym(end,2));

%% Animate side by side with energy plots
figure('Color','w','Position',[100 100 1200 500]);

subplot(2,2,1); hold on; axis equal; xlim([-2,2]); ylim([-2,2]);
xlabel('x'); ylabel('y'); title('Explicit Euler'); grid on
particle_exp = plot(X_exp(1,1), X_exp(1,2), 'ro','MarkerFaceColor','r','MarkerSize',8);
trail_exp = plot(X_exp(1,1), X_exp(1,2), 'b-','LineWidth',1);

subplot(2,2,2); hold on; axis equal; xlim([-2,2]); ylim([-2,2]);
xlabel('x'); ylabel('y'); title('Symplectic Euler'); grid on
particle_sym = plot(X_sym(1,1), X_sym(1,2), 'ro','MarkerFaceColor','r','MarkerSize',8);
trail_sym = plot(X_sym(1,1), X_sym(1,2), 'b-','LineWidth',1);

subplot(2,2,3); hold on; grid on
xlabel(' t'); ylabel(' $\mathcal{H}$ ','Interpreter','latex')

title('Explicit Euler Energy')
E_plot_exp = plot(t(1), E_exp(1), 'b-','LineWidth',1.5);

subplot(2,2,4); hold on; grid on
xlabel(' t'); ylabel(' $\mathcal{H}$ ','Interpreter','latex')

title('Symplectic Euler Energy')
E_plot_sym = plot(t(1), E_sym(1), 'b-','LineWidth',1.5);

for n = 1:5:length(t)
    % Update Explicit Euler
    subplot(2,2,1)
    set(particle_exp,'XData',X_exp(n,1),'YData',X_exp(n,2))
    set(trail_exp,'XData',X_exp(1:n,1),'YData',X_exp(1:n,2))
    
    subplot(2,2,3)
    set(E_plot_exp,'XData',t(1:n),'YData',E_exp(1:n))
    
    % Update Symplectic Euler
    subplot(2,2,2)
    set(particle_sym,'XData',X_sym(n,1),'YData',X_sym(n,2))
    set(trail_sym,'XData',X_sym(1:n,1),'YData',X_sym(1:n,2))
    
    subplot(2,2,4)
    set(E_plot_sym,'XData',t(1:n),'YData',E_sym(1:n))
    
    drawnow
end
