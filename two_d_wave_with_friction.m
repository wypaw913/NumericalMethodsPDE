%Math 511
%PDE HW 6, two-dimensional 2nd order wave equation with friction
%U_tt + eta*U_t = (c^2)(U_xx + U_yy)
%2nd order explicit FD scheme
%zero boundary conditions 
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

%wave speed, C
c = 1;

%coefficient of friction
eta = 0;

%resolution/spacing/dt
dt = 0.01;

%final time 0 for IC
Tf = 5;

%number of timesteps
timesteps = ceil(Tf/dt)+1;

%CFL (stability) constraint analogous to mu
%must be <=(1/2)*dx*(1/c)
r = c*dt/dx;

%constant involving friction and time for convenience

H = eta*dt/2;
B = (H-1)/(1+H);

%x vector of interior points from dx to b-dx with spacing dx
x_1=(dx:dx:b-dx);

%y vector of interior points from dy to d-dy with spacing dy
y_1=(dy:dy:d-dy);

%mesh grid for IC
[X_1,Y_1] = meshgrid(x_1,y_1);

%x vector of including BC with spacing dx
x =(a:dx:b);

%y vector including BC with spacing dy
y=x;

%mesh grid for surf
[X,Y] = meshgrid(x,y);

%%%%Building Matrix A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create center diagonal vector to use in A
center_diag = zeros( (N-1)*(M-1), 1) + 2 - 4*(r^2);

%create square matrix with center diagonal
A_center = diag(center_diag);

% create d+1 diagonal vector to use in A
diag_1 = zeros( (N-1)*(M-1)-1, 1);

for i=1 : length(diag_1)
    
    if mod(i, N-1) == 0
        
        diag_1(i) = 0;
        
    else
        diag_1(i) = r^2;
    end
end


% create d+N diagonal vector to use in A
diag_upper = zeros( (N-1)*(M-1)-(N-1), 1) + r^2 ;

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

%initial condition f(x,y) of interior pts only
%u_initial = 20.*X_1.*Y_1.*(X_1-1).*(Y_1-1);

%multiple regular waves
%u_initial = sin(X_1*10) + cos(Y_1*10);

%parallel ripples
%u_initial = cos(2*pi*(X_1+3*Y_1));

%two semicircles on boundary (double slit experiment)
%u_initial= 2*((X_1-.75).^2+(Y_1-0).^2 - .1^2 <0)+...
   %2*((X_1-.25).^2+(Y_1-0).^2 - .1^2 <0);
   
%quantum wave function
u_initial = sin(4*pi.*X_1).*sin(4*pi.*Y_1);

%jagged wall
%u_initial =  1.2^(-X_1.^2 - Y_1.^2);

%make a column vector
u_initial=u_initial(:);

%initial velocity g(x) interior pts
init_vel = 0*X_1.*Y_1;
%make a column vector
init_vel = init_vel(:);

%%%%Solve Linear System of Equations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
for i=0 : timesteps
        
    %preserve IC    
    if i==0
        
    u_new = u_initial;
    
    u_current = u_initial;
    
    %handle first timelevel
    elseif i==1
      
         %second order method: U1 = (1/2)A*U0 + ∆tg(x)
         u_1 = (1/(1-B))*((1/(1+H))*A*u_initial - ...
             (2*B*dt*init_vel) );
            
         %reassign
         u_current = u_initial;
         
    elseif i==2
        
        %Solving the linear system Un+1 = AUo - Un-1
        u_new = (A*u_1 +(H-1)*u_last)*(1/(H+1));
            
    %compute next time level
    else
        
        %Solving the linear system Un+1 = AUo - Un-1
        u_new = (A*u_current +(H-1)*u_last)*(1/(H+1));
    
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
    
%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    
    mesh(X, Y, u_complete)
    xlabel('x')
    ylabel('y')
    zlabel('u(x,y,t)')
    zlim([-2, 2])
    title(['Explicit Sol 2D Wave Eq. , r='...
        ,num2str(r), ', time Tf= ',num2str((i-1) * dt)...
        , ', J= ', num2str(N), ', Exec: ', num2str(toc)])
    
    
    %update U's
    u_last = u_current;
    
    u_current = u_new;

end

%figure(2)
%contour(X, Y, u_complete)


toc;