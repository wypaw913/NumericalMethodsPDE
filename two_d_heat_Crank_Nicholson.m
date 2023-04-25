%Math 511
%PDE HW 4, two-dimensional heat equation U_t = U_xx + U_yy
%2nd order Crank-Nicholson FD scheme, central in space in x and y
%zero and non-zero boundary conditions 
%compatible with any square mesh size on a square domain
clc; clear; close all;

tic;

% x domain [a, b]
a = 0;
b = 1;

% y domain [c, d]
c = 0;
d = 1;

%subintervals in x
N = 50;

%subintervals in y
M = N;

%resolution/spacing/dx
dx=(b-a)/N;

%resolution/spacing/dy
dy=(d-c)/M; 

%resolution/spacing/dt
dt = 0.1;

%final time 0 for IC
Tf = 0;

%number of timesteps
timesteps = Tf/dt;

%constant mu in x    
Vx = dt/(dx^2);

%constant mu in y
Vy = dt/(dy^2);

%define mu for square mesh
mu = Vx;

%x vector of interior points from dx to b-dx with spacing dx
x_1=(dx:dx:b-dx);

%y vector of interior points from dy to d-dy with spacing dy
y_1=(dy:dy:d-dy);

%mesh grid for surf
[X_1,Y_1] = meshgrid(x_1,y_1);

%x vector of interior points from dx to b-dx with spacing dx
x=(a:dx:b);

%y vector of interior points from dy to d-dy with spacing dy
y=(c:dy:d);

%mesh grid for surf
[X,Y] = meshgrid(x,y);

%%%%Building Matrix B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create center diagonal vector to use in B
center_diag = zeros( (N-1)*(M-1), 1) + 1-2*mu;

%create square matrix with center diagonal
B_center = diag(center_diag);

% create d+1 diagonal vector to use in B
diag_1 = zeros( (N-1)*(M-1)-1, 1);

for i=1 : length(diag_1)
    
    if mod(i, N-1) == 0
        
        diag_1(i) = 0;
        
    else
        diag_1(i) = mu/2;
    end
end


% create d+N diagonal vector to use in B
diag_upper = ones( (N-1)*(M-1)-(N-1), 1) * mu/2 ;

%create square matrix with d+N diagonal
B_3 = diag(diag_upper, N-1);

%transpose it to get square matrix with d-N diagonal
B_4 = diag(flip(diag_upper), -(N-1));

%create square matrix with d+1 diagonal
B_1 = diag(diag_1, 1);

%transpose it to get square matrix with d-1 diagonal
B_2 = diag(flip(diag_1), -1);

%add the 2 matrices to get d+1, d-1, d+N, d-N in one matrix
B_diags = B_1 + B_2 + B_3 + B_4;

%add the center diagonal to complete the matrix
B = B_diags + B_center;
    
%%%%Matrix B is fully formed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Building Matrix A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create center diagonal vector to use in A
center_diag = zeros( (N-1)*(M-1), 1) + 1+2*mu;

%create square matrix with center diagonal
A_center = diag(center_diag);

% create d+1 diagonal vector to use in A
diag_1 = zeros( (N-1)*(M-1)-1, 1);

for i=1 : length(diag_1)
    
    if mod(i, N-1) == 0
        
        diag_1(i) = 0;
        
    else
        diag_1(i) = -mu/2;
    end
end


% create d+N diagonal vector to use in A
diag_upper = ones( (N-1)*(M-1)-(N-1), 1) * -mu/2 ;

%create square matrix with d+N diagonal
A_3 = diag(diag_upper, N-1);

%transpose it to get square matrix with d-N diagonal
A_4 = diag(flip(diag_upper), -(N-1));

%create square matrix with d+1 diagonal
A_1 = diag(diag_1, 1);

%transpose it to get square matrix with d-1 diagonal
A_2 = diag(flip(diag_1), -1);

%add the 2 matrices to get d+1, d-1, d+N, d-N in one matrix
A_diags = A_1 + A_2 + A_3 + A_4;

%add the center diagonal to complete the matrix
A = A_diags + A_center;
    
%%%%Matrix A is fully formed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Create Initial Condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initial condition of interior pts only
u_initial = 20*X_1.*Y_1.*(1-X_1).*(1-Y_1);

%IC cylinder radius 0.25, height 2
%u_initial= 2*((X_1-.5).^2+(Y_1-.5).^2 - .25^2 <0);

%make a column vector
u_initial=u_initial(:);

%%%%Solve Linear System of Equations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_old = u_initial;
    
for i=0 : timesteps
        
        
    if i==0
    
            u_new = u_initial;

    %preserve IC   
    else
        
        %compute B*Un
        f = B*u_old;
        %Solving the linear system AUn+1 = BUn using Thomas Algo
        %u_new = tridiagonal(A, f);
        u_new = A\f;
    
    end
    
%%%%Add in Boundary Conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %so u_new can be turned back to col vector
    u_temp = u_new;
    
    %convert U from a column vector back into a matrix
    u_new = reshape(u_new, M-1, N-1);

    %pad the matrix with a border of zeros
    u_complete = padarray(u_new, [1,1]);
    
    %restore u_new as column vector
    u_new = u_temp;
    
    %plot
    figure(1)
    subplot(1, 2, 1)
    mesh(X, Y, u_complete)
    %pcolor(X, Y, u_complete)
    xlabel('x')
    ylabel('y')
    zlabel('u(x,y,t)')
    %zlim([0, 2.5])
    title(['Crank-Nicholson Sol 2D Heat Eq. , dt = '...
        ,num2str(dt), ', time Tf= ',num2str(i * dt)...
        , ', J= ', num2str(N), ', Exec: ', num2str(toc)])
    
    subplot(1,2,2)
    pcolor(X, Y, u_complete)
    colorbar
    xlabel('x')
    ylabel('y')
    zlabel('u(x,y,t)')
    
    %update u_old
    u_old = u_new;
end

toc;
