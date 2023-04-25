%this function is designed to implement the Thomas Algorithm to
%solve the vector equation x = A\p where A is a tridiagonal banded matrix
%with diagonals +-d
%written by Wyatt Paul using as a reference:
%"Tridiagonal Matrices: Thomas Algorithm" 
%by W.T. Lee of the University of Limerick
%modified by Wellington Woodward to handle spaced bands

function x = thomas_Paul_Woodward(A,p,d)

    %finds matrix dimension    
    n = size(A, 1);
    %initialize vector x, our solution
    x = zeros(n,1);
    %get diagonals from matrix as vectors
    b_diag = diag(A);
    c_diag = diag(A,d);
    a_diag = diag(A, -d);
    
    %Stage #1 forward sweep
    for i=1 : n
        
        %first row
        if i<=d
            
            c_diag(i) = c_diag(i)/ b_diag(i);
            p(i) = p(i)/ b_diag(i);
            b_diag(i) = 1;
            
        
        %interior rows
        elseif d<i && i<=n-d
            
            c_diag(i) = c_diag(i)/...
                (b_diag(i)-a_diag(i-d)*c_diag(i-d));
            p(i) = (p(i)-a_diag(i-d)*p(i-d))/...
                (b_diag(i)-a_diag(i-d)*c_diag(i-d));
            b_diag(i) = 1;
            a_diag(i-d) = 0;
            
        
        %last row
        else
           
            p(i) = (p(i)-a_diag(i-d)*p(i-d))/...
                (b_diag(i)-a_diag(i-d)*c_diag(i-d));
            b_diag(i) = 1;
            a_diag(i-d) = 0;
        end   
    end
    
    %Stage 2 backward substitution
    for i=n : -1: 1
        
        %first row
        if i<=d
            
            x(i) = p(i)-c_diag(i)*x(i+d);
            
        
        %interior rows
        elseif d<i && i<=n-d
            
            x(i) = p(i)-c_diag(i)*x(i+d);
           
        
        %last row do nothing
        else
            
            x(i) = p(i);
            
        end   
    end
end
