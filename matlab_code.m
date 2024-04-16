% Parameters
Lx = 1;  % Length of the domain in x-direction
Ly = 1;  % Length of the domain in y-direction
Nx = 50; % Number of grid points in x-direction
Ny = 50; % Number of grid points in y-directionn
T = 1;   % Total time
Nt = 100; % Number of time steps
alpha = 0.01; % Thermal diffusivity

% Grid spacing
dx = Lx / Nx;
dy = Ly / Ny;
dt = T / Nt;

% Create grid
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

% Initial condition
u0 = sin(pi*X) .* sin(pi*Y);

% Boundary conditions
leftBC = zeros(Ny, 1); % u(0, y, t) = 0
rightBC = zeros(Ny, 1); % u(Lx, y, t) = 0
topBC = zeros(1, Nx); % u(x, Ly, t) = 0
bottomBC = zeros(1, Nx); % u(x, 0, t) = 0

% Initialize solution matrix
u = zeros(Ny, Nx, Nt+1);
u(:,:,1) = u0;

% Finite difference scheme
for n = 1:Nt
    for i = 2:Nx-1
        for j = 2:Ny-1
            u(j, i, n+1) = u(j, i, n) + alpha * dt * ((u(j, i+1, n) - 2*u(j, i, n) + u(j, i-1, n)) / dx^2 + (u(j+1, i, n) - 2*u(j, i, n) + u(j-1, i, n)) / dy^2);
        end
    end

    % Apply boundary conditions
    u(:,1,n+1) = leftBC;
    u(:,end,n+1) = rightBC;
    u(1,:,n+1) = bottomBC;
    u(end,:,n+1) = topBC;
end

% Plotting
figure;
for i = 1:Nt+1
    surf(x, y, u(:,:,i));
    title(['Temperature distribution at time t = ', num2str((i-1)*dt)]);
    xlabel('x');
    ylabel('y');
    zlabel('Temperature');
    axis([0 Lx 0 Ly -1 1]); % Adjust axis limits if necessary
    pause(0.1); % Pause to display animation
end
