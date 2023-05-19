
A1=[4 1 -1;2 7 1;1 -3 12]; % The coefficient matrix for (i) part 
B1=[3 ;19;31]; % The RHS matrix for (i) part 
A2=[1 2 3;2 -1 2;3 1 -2]; % The coefficient matrix for (ii) part 
B2=[5  ;1 ;-1]; % The RHS matrix for (ii) part 
Xi1 = [0,0,0] ; % The initial approximation for (i) part 
Xi2 = [0,0,0] ;% The initial approximation for (ii) part 
X1=gauss_jacobi(A1,B1,Xi1); % X1 is the answer matrix to the (i) part
X2=gauss_jacobi(A2,B2,Xi2); % X1 is the answer matrix to the (ii) part
fprintf("The answer to part (a) is:\n x1=%.4f \t x2=%.4f \t x3=%.4f \n",X1(1),X1(2),X1(3));
fprintf("The answer to part (b) is:\n x1=%.4f \t x2=%.4f \t x3=%.4f \n",X2(1),X2(2),X2(3));
fprintf("The Gauss-Jacobi method is guaranteed to converge for part (a) since the matrix in question is diagonally \n" + ...
    "dominant (sufficient condition) while it fails to converge in case of part(b) with any initial approximation \n" + ...
    "(unless we start with the solution of the system of equations as the initial value).\n")
function X2=gauss_jacobi(A,B,Xi)
    k=size(A); 
    n=k(1);
    X1=Xi ; % 1st approximation
    X2=Xi + 1 ;
    while(check(X2,X1,n)>=0.000001)
        X1=X2;
        for i=1:n
            temp=B(i);
            for j=1:n
                if(i~=j)
                    temp=temp-A(i,j)*X1(j); % Changing the value in the next iteration
                end
            end
            X2(i)=temp/A(i,i); 
        end  
    end    
end


function c=check(A,B,n)
% This function is used to define the terminating condition of the loop
% When the 2 consecutive iterations are too 'close' we stop the process.
    c=0;
    for i=1:n
        c=c+abs(A(i)-B(i));
    end
end

