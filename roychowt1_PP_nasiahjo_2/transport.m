function phi_new = transport(phi,u, N, dt, dx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
size_phi = size(phi);
size_u = size(u);

if size_phi(1,1) ~= 129
    error('Matrix of phi is not the correct size');
    if size_phi(1,2) ~= 1
        error('Matrix of phi is not the correct size');
    end
end

if size_u(1,1) ~= 129
    error('Matrix of phi is not the correct size');
    if size_u(1,2) ~= 1
        error('Matrix of phi is not the correct size');
    end
end

A = zeros(N+1);
A(1,1) = -1;
A(1,2) = 1;
A(N+1,N+1) =1;
A(N+1,N) = -1;
b = zeros(129,1);
for n = 2:N
    A(n,n-1) = (-dt)*u(n,1);
    A(n,n) = 1 + (dt*(u(n+1,1) - u(n-1,1)))/(2*dx);
    A(n,n+1) = dt*u(n,1)/(2*dx);
    b(n,1) = phi(n,1);
end

phi_new = A\b;   
end



