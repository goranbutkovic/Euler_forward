%% 2D particle in nonlinear potential (Hénon–Heiles): Explicit vs Symplectic Euler
clear; clc; close all;

% Parameters
m = 1; kx = 1; ky = 1; lambda = 0.3;
dt = 0.01; tf = 100; t = 0:dt:tf;

% Potential and derivatives
V = @(x,y) 0.5*kx*x.^2 + 0.5*ky*y.^2 + lambda*(x.^2.*y - (1/3)*y.^3);
dVdx = @(x,y) kx*x + 2*lambda*x.*y;
dVdy = @(x,y) ky*y + lambda*(x.^2 - y.^2);

% Initial conditions
x0 = -0.4; y0 = 3; vx0 = 0; vy0 = 0;

%% Explicit Euler
X_exp = zeros(length(t),4); X_exp(1,:) = [x0,y0,vx0,vy0]; E_exp = zeros(length(t),1);

for n = 1:length(t)-1
    x = X_exp(n,1); y = X_exp(n,2); vx = X_exp(n,3); vy = X_exp(n,4);
    Fx = -dVdx(x,y); Fy = -dVdy(x,y);
    vx_new = vx + dt*Fx/m; vy_new = vy + dt*Fy/m;
    x_new = x + dt*vx; y_new = y + dt*vy;
    X_exp(n+1,:) = [x_new, y_new, vx_new, vy_new];
    E_exp(n) = 0.5*m*(vx^2+vy^2) + V(x,y);
end
E_exp(end) = 0.5*m*(X_exp(end,3)^2+X_exp(end,4)^2)+V(X_exp(end,1),X_exp(end,2));

%% Symplectic Euler
X_sym = zeros(length(t),4); X_sym(1,:) = [x0,y0,vx0,vy0]; E_sym = zeros(length(t),1);

for n = 1:length(t)-1
    x = X_sym(n,1); y = X_sym(n,2); vx = X_sym(n,3); vy = X_sym(n,4);
    Fx = -dVdx(x,y); Fy = -dVdy(x,y);
    vx_new = vx + dt*Fx/m; vy_new = vy + dt*Fy/m;
    x_new = x + dt*vx_new; y_new = y + dt*vy_new;
    X_sym(n+1,:) = [x_new, y_new, vx_new, vy_new];
    E_sym(n) = 0.5*m*(vx_new^2 + vy_new^2) + V(x_new,y_new);
end
E_sym(end) = 0.5*m*(X_sym(end,3)^2+X_sym(end,4)^2)+V(X_sym(end,1),X_sym(end,2));

%% Animate Explicit Euler (fixed view)
fig1 = figure('Color','w','Position',[100 100 800 600]); 
hold on; axis equal;

% Fixed axis limits
xlim([-2.7 2.7]); 
ylim([-2.5 3]); 

grid on;
xlabel('x'); ylabel('y');

% Initial particle and trail
particle_exp = plot(X_exp(1,1), X_exp(1,2), 'ro','MarkerFaceColor','r','MarkerSize',8);
trail_exp = plot(X_exp(1,1), X_exp(1,2), 'b-','LineWidth',1);

% Animation loop
for n = 1:5:length(t)
    set(particle_exp,'XData',X_exp(n,1),'YData',X_exp(n,2));
    set(trail_exp,'XData',X_exp(1:n,1),'YData',X_exp(1:n,2));
    drawnow; pause(0.01);
end

%% Animate Symplectic Euler (same fixed view)
fig2 = figure('Color','w','Position',[150 150 800 600]); 
hold on; axis equal;

% Fixed axis limits same as Explicit Euler
xlim([-2.7 2.7]); 
ylim([-2.5 3]); 

grid on;
xlabel('x'); ylabel('y');

particle_sym = plot(X_sym(1,1), X_sym(1,2), 'ro','MarkerFaceColor','r','MarkerSize',8);
trail_sym = plot(X_sym(1,1), X_sym(1,2), 'b-','LineWidth',1);

for n = 1:5:length(t)
    set(particle_sym,'XData',X_sym(n,1),'YData',X_sym(n,2));
    set(trail_sym,'XData',X_sym(1:n,1),'YData',X_sym(1:n,2));
    drawnow; pause(0.01);
end


%% Energy plot Explicit Euler
fig3 = figure('Color','w','Position',[200 200 800 600]);
plot(t,E_exp,'b','LineWidth',1.5); grid on;
xlabel('t'); ylabel(' $\mathcal{H}$ ','Interpreter','latex'); %title('Energy - Explicit Euler');

%% Energy plot Symplectic Euler
fig4 = figure('Color','w','Position',[250 250 800 600]);
plot(t,E_sym,'b','LineWidth',1.5); grid on;
xlabel('t'); ylabel(' $\mathcal{H}$ ','Interpreter','latex'); %title('Energy - Symplectic Euler');
