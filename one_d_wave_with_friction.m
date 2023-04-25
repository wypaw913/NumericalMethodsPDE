%Math 511
%PDE HW 6, one-dimensional 2nd order wave equation with friction
%U_tt + eta*U_t = (c^2)U_xx
%2nd order explicit FD scheme
%zero boundary conditions 
clc; clear; close all;

tic;

% x domain [a, b]
a = 0;
b = 1;

%subintervals in x
N = 100;

%resolution/spacing/dx
dx=(b-a)/N;

%x vector of all points with spacing dx
x=(a:dx:b);

%x vector of interior points from dx to b-dx with spacing dx
x_1=(dx:dx:b-dx);

%resolution/spacing/dt
dt = 0.01;

%final time 0 for IC
Tf = 3;

%wave speed, C
c = 1;

%coefficient of friction
eta = 0;

%CFL (stability) constraint analogous to mu
%must be <=1
r = c*dt/dx;

%constant involving friction and time for convenience
B = eta*dt/2;

%number of timesteps
timesteps = ceil(Tf/dt);


%%%%Building Matrix A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create center diagonal vector to use in A
center_diag = zeros( (N-1), 1) + (2*(1-r^2));

%create square matrix with center diagonal
A_center = diag(center_diag);

% create d+1 diagonal vector to use in A
diag_1 = zeros( N-2, 1) + r^2;

%create square matrix with d+1 diagonal
A_1 = diag(diag_1, 1);

%transpose it to get square matrix with d-1 diagonal
A_2 = diag(flip(diag_1), -1);

%add the 2 matrices to get d+1, d-1 in one matrix
A_diags = A_1 + A_2;

%add the center diagonal to complete the matrix
A = A_diags + A_center;
    
%%%%Matrix A is fully formed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Create Initial Conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initial condition of interior pts only

%initial displacement f(x)interior pts sine wave
u_initial = sin(2.*pi.*x_1./b);

%spike
%u_initial = exp(-400*(x_1-0.5).^2);

%make a column vector
u_initial = u_initial(:);


%initial velocity g(x) interior pts
init_vel = 0*x_1;
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
         u_1 = 0.5*A*u_initial + dt*init_vel;
            
         %reassign
         u_current = u_initial;
         
         
    %compute next time level
    elseif i==2
      
         %second order method: U1 = (1/2)A*U0 + ∆tg(x)
         u_new = (A*u_1 +(B-1)*u_last)*(1/(B+1));
            
         %reassign
         u_current = u_initial;
         
         
    %compute next time level
    else
        
        %Solving the linear system Un+1 = AUo - Un-1
        u_new = (A*u_current +(B-1)*u_last)*(1/(B+1));
    
    end
    
%%%%Add in Boundary Conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    %convert u from col vect to horiz vector and add BC's to either end
    u_complete = [0, u_new', 0];
    
%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1)
    plot(x, u_complete, '-o')
    xlabel('x')
    ylabel('U(x,t)')
    ylim([-1,1])
    title(['Explicit Sol 1D Wave Eq. , r='...
        ,num2str(r), ', time Tf= ',num2str(i * dt)...
        , ', J= ', num2str(N), ', Exec: ', num2str(toc)])
    
    
    %update U's
    u_last = u_current;
    
    u_current = u_new;

end

toc;