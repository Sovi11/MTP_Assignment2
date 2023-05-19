
A=randn(10,10) ; 
b=randn(10,1) ;
a1 =solve(A,b) ;
% a1 is the matrix obtained using Doolittle decomposition
a2  = A\b; 
% a2 is the matrix obtained using inbuilt matrix 'division'
a = [a1 a2] ;
disp(" ")
disp(" The 10x1     The 10x1  ")
disp(" matrix by    matrix by")
disp(" Doolittle    A\b method")
disp("Decomposition")

disp(a) ;
function x = solve(A,b)
n = (size(A)) ; 
n= n(1) ;
Y=zeros(n,1) ;
rr =break_in_lu(A) ; % rr is an augmented matrix of L and U 
% breaking A into L*U where L is a lower triangular matrix and U is an
% upper triangular matrix
L = rr(:,1:n) ;
U = rr(:,n+1:2*n) ;
% Solve LY = b then 
% solve UX = Y hence the given X satisfies LUX = b hence AX = b
for i = 1:n
    sm = b(i) ;
    for j=1:i-1
        sm = sm - L(i,j)*Y(j) ;
    end 
    Y(i) = sm/L(i,i) ;
end 
x = zeros(n,1) ;
for i= n:-1:1
    sm = Y(i) ;
    for j= i+1:n
        sm =sm - U(i,j)*x(j) ;
    end
    x(i) = sm/U(i,i) ;
end 
end

function f = break_in_lu(A)
% Breaking the matrix into Lower and Upper triangular matrices L and U such
% that LU = A
n = (size(A)) ; 
n= n(1) ;
U = zeros([n,n]) ;
L = zeros([n,n]) ;
for i = 1:n
    for k= i:n 
        sm = 0 ; 
    for j = 1:i
        sm = sm + L(i,j)*U(j,k) ;
    end
    U(i,k) = A(i,k) - sm ;
    end 
    for k = i:n
        if (i ==k)
            L(i,i) = 1 ;    
        else 
            sm = 0 ; 
            for j = 1:i
                sm = sm + L(k,j)*U(j,i) ;
            end
            L(k,i) = (A(k,i) - sm)/U(i,i);
        end 
    end
end 
f=[L,U] ;
end 

 
