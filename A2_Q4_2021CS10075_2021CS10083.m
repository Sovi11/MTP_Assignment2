
data = [12, 279.2;8, 177.2;5,106.8  ] ; % Data table for time vs velocity
 A = zeros(3,3) ; % This is the coefficient matrix A in the equation AX = b
 for i = 1:3 
     A(i,1) = data(i,1)^2 ;
     A(i,2) = data(i,1) ;
     A(i,3) = 1 ;
 end 
B = data(:,2) ;% This is the RHS matrix b in the equation AX = b
a =1 ; 
b = 2 ;
c = 5; 
% X is [a ; b ; c] in the equation AX = b 
% This is gauss seidel method in which we change the variables in each
% iteration unlike the Guass jacobi method in which to calculate the k+1 th
% iteration kth iteration values for all variables are used
while (abs((a*A(1,1) + b*A(1,2) + c*A(1,3) - B(1))) > 10^(-6)) 
    % Terminating condition is when the equations are almost satisfied 
    a = (B(1) - A(1,2)*b - A(1,3)*c)/A(1,1) ;
    b = (B(2) - A(2,1)*a - A(2,3)*c)/A(2,2) ;
    c = (B(3) - A(3,2)*b - A(3,1)*a)/A(3,3) ;
end
fprintf("a = %.6f \nb = %.6f \nc = %.6f \n",a,b,c)

