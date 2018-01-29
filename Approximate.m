function [P] = Approximate(a, b, lambda, d, N)

%Creating our discrete domain as matrix P
%   As we are using matrices for a triangular domain there will be some
%   waste
P = zeros(N,N);

% handy constant values to keep expressions short in the following 
% approximations
h = (b-a)/N;
v = -(lambda + d);

%Applying the boundary conditions
%   Setting the right column to 0
P(:,N)=0;
%   Applying the diagonal condition
for i = 1:N
    x = a + h*i;
    P(i,i)= (-1/2)*v*(x - b);
end

% Going through the whole domain to compute the next approximation for each
% cell of the matrix. We do not touch the areas whose boundary conditions 
% was set.
for j = N-2:-1:1
    for i = N-1:-1:j+1
        if i == j+1
            % Points directly below the diagonal
            P(j,i) = (P(j,i+1) + P(j,i-1))/2;
        else
            % Regular approximation formula
            P(j,i) = (P(j+1,i+1) + P(j+1,i-1) - P(j+2,i))-v*h*h*P(j+1,i);
        end
    end
end

end

