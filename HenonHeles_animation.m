%% 2D čestica u kaotičnom Hénon–Heiles potencijalu s prikazom površine

clear all; close all; clc;

%% Parameters
m = 1;              % mass
kx = 1; ky = 1;     % quadratic spring constants
lambda = 0.3;       % nonlinearity strength (try 0.1–1.0)
dt = 0.01;          % time step
tf = 1000;            % final time
nsteps = tf/dt;

%% Potential function (Hénon–Heiles)
V = @(x,y) 0.5*kx*x.^2 + 0.5*ky*y.^2 + lambda*(x.^2.*y - (1/3)*y.^3);

%% Forces (negative gradient of potential)
Fx = @(x,y) - (kx*x + 2*lambda*x.*y);
Fy = @(x,y) - (ky*y + lambda*(x.^2 - y.^2));

%% Initial conditions
x = zeros(nsteps,1);
y = zeros(nsteps,1);
vx = zeros(nsteps,1);
vy = zeros(nsteps,1);

x(1) = -0.4;     % initial x
y(1) = 3;     % initial y
vx(1) = 0;   % initial velocity x
vy(1) = 0;    % initial velocity y

%% Integrate using symplectic Euler
for i = 1:nsteps-1
    % Update velocities (p_{n+1} = p_n + dt * F(x_n, y_n))
    vx(i+1) = vx(i) + dt/m * Fx(x(i),y(i));
    vy(i+1) = vy(i) + dt/m * Fy(x(i),y(i));
    
    % Update positions (q_{n+1} = q_n + dt * v_{n+1})
    x(i+1) = x(i) + dt * vx(i+1);
    y(i+1) = y(i) + dt * vy(i+1);
end

%% Compute potential height for visualization
z = V(x,y);

%% Plot potential surface with exact particle path
[Xgrid,Ygrid] = meshgrid(-5:0.03:5, -5:0.03:5);
Vgrid = V(Xgrid,Ygrid);

figure('Color','w')
surf(Xgrid,Ygrid,Vgrid,'FaceAlpha',0.6,'EdgeColor','none')
colormap(parula)
hold on

%% Particle animation with growing trajectory
particle = plot3(x(1),y(1),z(1),'ro','MarkerFaceColor','r','MarkerSize',8);
traj = plot3(x(1),y(1),z(1),'r','LineWidth',1.5);  % start with first point

xlabel('x'); ylabel('y'); zlabel('V(x,y)');
%title('Particle moving in chaotic Hénon–Heiles potential')
view(45,30)
grid on

for i = 2:5:nsteps
    % update particle position
    set(particle,'XData',x(i),'YData',y(i),'ZData',z(i))
    
    % update trajectory
    set(traj,'XData',x(1:i),'YData',y(1:i),'ZData',z(1:i))
    
    drawnow
end
