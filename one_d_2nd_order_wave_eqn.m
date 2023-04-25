%Math 511
%PDE HW 5, one-dimensional 2nd order wave equation U_tt = (c^2)U_xx
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
dt = 0.008;

%final time 0 for IC
Tf = 1;

%wave speed, C
c = 1;

%CFL (stability) constraint analogous to mu
%must be <=1
r = c*dt/dx;

%number of timesteps
timesteps = ceil(Tf/dt);

%%%%for u_reference%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%resolution/spacing/dt
dt_ref = 1e-5;

%CFL (stability) constraint analagous to mu
%must be <=1
r_ref = c*dt_ref/dx;

%number of timesteps
timesteps_ref = ceil(Tf/dt_ref);

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

%%%%Building Matrix A_ref%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create center diagonal vector to use in A
center_diag = zeros( (N-1), 1) + (2*(1-r_ref^2));

%create square matrix with center diagonal
A_center = diag(center_diag);

% create d+1 diagonal vector to use in A
diag_1 = zeros( N-2, 1) + r_ref^2;

%create square matrix with d+1 diagonal
A_1 = diag(diag_1, 1);

%transpose it to get square matrix with d-1 diagonal
A_2 = diag(flip(diag_1), -1);

%add the 2 matrices to get d+1, d-1 in one matrix
A_diags = A_1 + A_2;

%add the center diagonal to complete the matrix
A_ref = A_diags + A_center;
    
%%%%Matrix A_ref is fully formed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Create Initial Condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initial condition of interior pts only

%initial displacement f(x)interior pts spike
%u_initial = exp(-400*(x_1-0.5).^2);
%make a column vector
%u_initial = u_initial(:);

%initial displacement f(x)interior pts square wave
%u_initial = (0.4<=x_1)&(x_1<=0.6);
%make a column vector
%u_initial = u_initial(:);

%initial displacement f(x)interior pts sine wave
u_initial = sin(pi*x_1);
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
       
         %first order method: U1 = U0 + ∆tg(x)
         %u_new = u_initial + dt*init_vel;


         %second order method: U1 = (1/2)A*U0 + ∆tg(x)
         u_new = 0.5*A*u_initial + dt*init_vel;
            
         %reassign
         u_current = u_initial;
            
    %compute next time level
    else
        
        %Solving the linear system Un+1 = AUo - Un-1
        u_new = A*u_current - u_last;
    
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

% %%%%compute u reference%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=0 : timesteps_ref
%         
%     %preserve IC 
%     if i==0
%        
%         u_new_ref = u_initial;
%     
%         u_current_ref = u_initial;
%         
%     %%handle first time level   
%     elseif i==1   
%     
%         %first order method: U1 = U0 + ∆tg(x)
%         %u_new_ref = u_initial + dt_ref*init_vel;
% 
% 
%         %second order method: U1 = (1/2)A*U0 + ∆tg(x)
%         u_new_ref = 0.5*A_ref*u_initial + dt_ref*init_vel;
% 
%         %reassign U's
%         u_current_ref = u_initial;
%     
%     %compute next time level
%     else
%         
%         %Same for u reference
%         u_new_ref = A_ref*u_current_ref - u_last_ref;
%     
%     end
%     
% %%%%Add in Boundary Conditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% 
%     %same for u reference
%     u_ref = [0, u_new_ref', 0];
%       
%     %reassign u's
%     %same for u reference
%     u_last_ref = u_current_ref;
%     
%     %same for u reference
%     u_current_ref = u_new_ref;
% end

%%compute max solution error at Tf%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%max_error = max(abs(u_ref - u_complete))

toc;