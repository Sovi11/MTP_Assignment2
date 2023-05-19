
A=[9 3 2 0 7;7 6 9 6 4;2 7 7 8 2;0 9 7 2 2;7 3 6 4 3]; % coefficient matrix
b=[35; 58; 53; 37; 39]; % RHS matrix 
X=gauss_elimination(A,b); %Solution vector X 
% To solve AX =b 
fprintf("x1 = %.4f \nx2 = %.4f \nx3 = %.4f \nx4 = %.4f \nx5 = %.4f\n",X(1),X(2),X(3),X(4),X(5));
 
function x = gauss_elimination(A,b)
    A=[A b]; %Augmented Matrix
    [n, p]= size(A); %Dimensions of Augmented Matrix
    
    for k = 1:n-1
        [index,hm]=max(A(k:n,k)); %Finding the index of the maximum value among the values in the column below A(k,k)
        A([k hm+k-1],:)=A([hm+k-1 k],:); %Pivoting (swapping the rows of the augmented matrix)
   
            
        %Forward Elimination
        for i = k+1:n
            m = A(i,k) / A(k,k);
            for j=k:p
                A(i,j)=A(i,j)-m*A(k,j); % elimination , making the coefficient matrix upper triangular
            end
         end
    end
    

    % Backward Substitution (solving for the values of X)
    x = zeros(n,1);
    x(n) = A(n,p) / A(n,n);
    for i = n-1:-1:1
        x(i) = (A(i,p) - A(i,i+1:p-1)*x(i+1:p-1)) / A(i,i); %Calculating every x[i] using back substitution
    end
end

