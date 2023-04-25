%Math 511
%PDE hw 3, Laplace's Equation Uxx + Uyy =0
%steady-state behavior of 2d heat equation
%2nd order FD scheme, central in space in x and y
%zero and non-zero boundary conditions 
%compatible with any square mesh size on a square domain
clc; clear; close all;

% x domain [a, b]
a = 0;
b = 0.5;

% y domain [c, d]
c = 0;
d = 0.5;

%subintervals in x
N = 100;

%subintervals in y
M = N;

%resolution/spacing/dx
dx=(b-a)/N;

%resolution/spacing/dy
dy=(d-c)/M;            

%x vector from a to b with spacing dx
x=(a:dx:b);

%y vector from c to d with spacing dy
y=(c:dy:d);

%mesh grid for surf
[X,Y] = meshgrid(x,y);


%%%%Building Matrix A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create center diagonal vector to use in A
center_diag = zeros( (N-1)*(M-1), 1) -4;

%create square matrix with center diagonal
A_center = diag(center_diag);

% create d+1 diagonal vector to use in A
diag_1 = zeros( (N-1)*(M-1)-1, 1);

for i=1 : length(diag_1)
    
    if mod(i, N-1) == 0
        
        diag_1(i) = 0;
        
    else
        diag_1(i) = 1;
    end
end


% create d+N diagonal vector to use in A
diag_upper = ones( (N-1)*(M-1)-(N-1), 1);

%create square matrix with d+3 diagonal
A_3 = diag(diag_upper, N-1);

%transpose it to get square matrix with d-3 diagonal
A_4 = diag(flip(diag_upper), -(N-1));

%create square matrix with d+1 diagonal
A_1 = diag(diag_1, 1);

%transpose it to get square matrix with d-1 diagonal
A_2 = diag(flip(diag_1), -1);

%add the 2 matrices to get d+1, d-1, d+3, d-3 in one matrix
A_diags = A_1 + A_2 + A_3 + A_4;

%add the center diagonal to complete the matrix
A = A_diags + A_center;
    
%%%%Matrix A is fully formed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Building Vector b%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create Left BC vector
left_BC = zeros( (N-1)*(M-1), 1);

%create reduced M+1 vector to use in solution
red_left_BC = zeros(M+1, 1);

%variable to increment algorithm
n=0;

for i=1 : length(left_BC)
    
    %using an algorithm to generate L based on N subintervals in x
    if mod(i, (n*(N-1))+1 ) == 0
        
        %Multiply BC by y vector element
        left_BC(i) = (-1*y(n+1)) * 0;
        %add to reduced vector
        red_left_BC(n+1) = abs(left_BC(i));
        
        %increment algo counter
        n = n+1;        
    end
end

% create Right BC vector
right_BC = zeros( (N-1)*(M-1), 1);

%create reduced M+1 vector to use in solution
red_right_BC = zeros(M+1, 1);

%variable to increment algorithm
n=1;

for i=1 : length(right_BC)
    
    %using an algorithm to generate R based on N subintervals in x
    if (mod(i, n*(N-1)) == 0) && i~=0
        
        %multiply BC by y vector element
        right_BC(i) = (-1*y(n+1)) * 200;
        %add to reduced vector
        red_right_BC(n+1) = abs(right_BC(i));
        
        %increment algo counter
        n = n+1;        
    end
end

% create Bottom BC vector
bottom_BC = zeros( (N-1)*(M-1), 1);

%create reduced N+1 vector to use in solution
red_bot_BC = zeros(N+1, 1);

for i=1 : N-1
        
        %multiply BC by x vector element
        bottom_BC(i) = (-1* x(i+1)) * 0;
        %add to reduced vector
        red_bot_BC(n+1) = abs(bottom_BC(i));
              
end

% create Top BC vector
top_BC = zeros( (N-1)*(M-1), 1);

%create reduced N+1 vector to use in solution
red_top_BC = zeros(N+1, 1);

%to increment x vector index
n=1;

for i=(N-1)*(M-1)-(N-2) : (N-1)*(M-1)
        
        %multiply BC by x vector element
        top_BC(i) = (-1* x(n+1)) * 200;
        %add to reduced vector
        red_top_BC(n+1) = abs(top_BC(i));
        %increment algo counter
        n = n+1;       
end

%create vector b by summing all four bc vectors L,R,B,T
f = left_BC + right_BC + bottom_BC + top_BC;

%%%%Vector b is fully formed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Solve the linear system Au=f by u = A\f;
u = A\f;

%convert U from a column vector back into a matrix
u = reshape(u, M-1, N-1);

%pad the matrix with a border of zeros
u = padarray(u, [1,1]);

%create (M+1)x(N+1) matrix to house BC's
u_complete = zeros(M+1, N+1);

%make left row left BC
u_complete(:, 1) = red_left_BC;

%make right row right BC
u_complete(:, N+1) = red_right_BC;

%make bottom row bottom BC
u_complete(1, :) = red_bot_BC;

%make top row top BC
u_complete(N+1, :) = red_top_BC;

% add the cherry on top, top right corner BC
u_complete(M+1, N+1) = 200*x(length(x));

%add u to u_complete to finish u_complete
u_complete = u_complete + u;

%exact solution
%u_exact = 400.*X_fine.*Y_fine;
u_exact = 400.*X.*Y;

%%%%solution error%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_error = abs(u_exact - u_complete);

%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(2,2,1)
%numerical solution
surf(x,y, u_complete)
shading interp
colorbar;
xlim([0, 0.5])
ylim([0, 0.5])
xlabel('x')
ylabel('y')
zlabel('u(x,y,t)')
title(['Steady-State Soln to 2-D Heat Eqn with'...
' Non-Zero BCs, N=50'])

subplot(2,2,3)
pcolor(X, Y, u_complete)
shading interp
colorbar
xlabel('x')
ylabel('y')
title('N=50')

subplot(2,2,2)

%numerical solution
surf(X,Y, u_exact)
shading interp
colorbar;
xlim([0, 0.5])
ylim([0, 0.5])
xlabel('x')
ylabel('y')
zlabel('u(x,y,t)')
title(['Steady-State Soln to 2-D Heat Eqn with'...
' Non-Zero BCs, Exact Solution'])

subplot(2,2,4)
pcolor(X, Y, u_exact)
shading interp
colorbar
xlabel('x')
ylabel('y')
title('Exact')

%%solution error plotting
figure(2)
surf(X,Y, u_error)
shading interp
colorbar;
xlim([0, 0.5])
ylim([0, 0.5])
xlabel('x')
ylabel('y')
zlabel('u(x,y,t)')
title(['Solution Error, N=40'])