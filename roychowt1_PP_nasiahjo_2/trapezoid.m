function [lhs_integral] = trapezoid(f_n,N,dx)

%assigns the size of the each input variables an array
size_N=size(N);
size_f_n=size(f_n);
size_dx=size(dx);

%determines if the size of the array is equal to the specified value.
%if it isn't equal to the specified value the function returns error
if(size_N~=[1 1])
    error(error)
end
if(size_f_n~=[N+1 1])
    error(error)
end
if(size_dx~=[1 1])
    error(error)
end

%initializes the variable n and the while loop will loop N times and create
%a summation of the previous value of f_n and the current value then
%multiply it by dx
for n=1:1:N
    lhs_integral=((f_n(n,1))+(f_n(n+1,1)))/2*dx;
end
end

