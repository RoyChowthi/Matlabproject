%scalars
L=18;
g=9.81;
H=0.75;
rho=1000;
w=1.5;
N=128;
M=512;
dx=L/N;
dt=(0.9*dx)/sqrt(g*H);
%vectors and 1D arrays
x=linspace(0,N*dx,N+1);
theta=2*pi*x/L;
t=linspace(0,dt*M,M);
%matrices and 2D arrays
u=zeros(N+1,M);
h=zeros(N+1,M);
eta=zeros(N+1,M);
phi=zeros(N+1,M);
u(:,1)=(5/11)*sqrt(g*H)*sin(theta/2);
h(:,1)=(H)*((1/3).*cos((theta-pi)/2 ).^16 + (1/5).*(cos((theta-pi/2)/4).^32)) + (1-(x/L))+(691/1818);
eta(:,1)=h(:,1)-H;
phi(:,1)=(1/100)*exp(-128*((x./L)-(3/4)).^2);
%part a
%load
load h_lin.dat; 
load u_lin.dat;
load soln.mat;
figure(1)
subplot(3,1,1);
plot(x,u_lin(:,1))
xlabel('x')
ylabel('u lin')
title('u lin vs x')
subplot(3,1,2);
plot(x,h_lin(:,1))
xlabel('x')
ylabel('h lin')
title('h lin vs x')
subplot(3,1,3);
plot(x,phi(:,1))
xlabel('x')
ylabel('phi')
title('phi vs x')
%part b
figure(2)
subplot(3,1,1);
plot(x,u_nlin(:,1))
xlabel('x')
ylabel('u nlin')
title('u nlin vs x')
subplot(3,1,2);
plot(x,h_nlin(:,1))
title('h lin vs x')
xlabel('x')
ylabel('h lin')
subplot(3,1,3);
plot(x,phi_nlin(:,1))
title('phi nlin vs x')
xlabel('x')
ylabel('phi nlin')
%part c
figure(3)
subplot(5,2,1);
plot(x,h_lin(:,1))
title('m=1')
xlabel('x');
ylabel('h');
subplot(5,2,2)
plot(x,h_nlin(:,1))
title('m=1')
xlabel('x');
ylabel('h');
subplot(5,2,3)
plot(x,h_lin(:,128))
title('m=128')
xlabel('x');
ylabel('h');
subplot(5,2,4)
plot(x,h_nlin(:,128))
title('m=128')
xlabel('x');
ylabel('h');
subplot(5,2,5)
plot(x,h_lin(:,256))
title('m=256')
xlabel('x');
ylabel('h');
subplot(5,2,6)
plot(x,h_nlin(:,256))
title('m=256')
xlabel('x');
ylabel('h');
subplot(5,2,7)
plot(x,h_lin(:,384))
title('m=384')
xlabel('x');
ylabel('h');
subplot(5,2,8)
plot(x,h_nlin(:,384))
title('m=384')
xlabel('x');
ylabel('h');
subplot(5,2,9)
plot(x,h_lin(:,512))
title('m=512')
xlabel('x');
ylabel('h');
subplot(5,2,10)
plot(x,h_nlin(:,512))
title('m=512')
xlabel('x');
ylabel('h');

%----------------------------------------------------------------------------------
%------------------------------Part 2 --------------------------
flag = 1;
while flag == 1
    if flag == 1
    selection = input('Please enter 0 for linear and 1 for non-linear . If either number is not 0 or one, you will be asked to do so. ');
        if selection == 1 
            flag = 0;
        end 
        if selection == 0 
            flag = 0;
        end 
        
    end     
end 

for m = 2:512
    if selection == 1
        %Why are the variables not updating correctly is the variables 
        
        [u(:,m),h(:,m),eta(:,m),phi(:,m)] = nonlinear(u(:,m-1),h(:,m-1),eta(:,m-1),phi(:,m-1),N,dx,dt,g,H);
        
        %each loop fill out data table 
%         u(:,i)= u_est(:,i);
%         h(:,i)= h_est(:,i);
%         eta(:,i)= eta_est(:,i);
%         phi(:,i)= phi_est(:,i);
    end 
    if selection == 0
       [u(:,m),h(:,m),eta(:,m)] = gravity(u(:,m-1),h(:,m-1),eta(:,m-1),N,dx,dt,g,H);
       phi(:,m) = transport(phi(:,m-1),u(:,m-1), N, dt, dx);
    end
end
%to plot both subplots run the program twice the first time type in 0 for
%linear and the second time type in 1 for nonlinear
if selection == 0
    figure(4)
    hold on
    subplot(1,2,1)
    plot(x,phi(:,1),'-b',x,phi(:,128),'--r',x,phi(:,184),':k');
    title('Linear Tracer Concentration');  
    xlabel('Length');
    ylabel('Tracer Concentration');
    legend('At m = 1','At m = 128','At m = 184');
end 

 if selection == 1
     figure(4)
     hold on
     subplot(1,2,2)
     plot(x,phi(:,1),'-b',x,phi(:,128),'--r',x,phi(:,184),':k');
     title('Nonlinear Tracer Concentration');  
     xlabel('Length');
     ylabel('Tracer Concentration');
     legend('At m = 1','At m = 128','At m = 184');
     
 end 
 hold off
%-----------------------------------Part 3----------------------------------

f_n=zeros(N+1,1);
Available_Energy=zeros(N+1,M);
Fluid_Mass=zeros(N+1,M);
Kinetic_Energy=zeros(N+1,M);
Potential_Energy=zeros(N+1,M);
for i=1:1:N+1
    for j=1:1:M
    %available Energy
    f_n(i,1)=((1/2)*h(i,j).*(u(i,j).^2))+((1/2)*g*(eta(i,j).^2));
    lhs_integral = trapezoid(f_n,N,dx);
    Available_Energy(i,j)=(rho*w)*lhs_integral;
    Available_Energy(i,j)=Available_Energy(i,j)/Available_Energy(i,1);
    end
end
for i=1:1:N+1
    for j=1:1:M
    %Fluid Mass
    f_n(i,1)=h(i,j);
    lhs_integral = trapezoid(f_n,N,dx);
    Fluid_Mass(i,j)=(rho*w)*lhs_integral;
    Fluid_Mass(i,j)=Fluid_Mass(i,j)/Fluid_Mass(i,1);
    end    
end
for i=1:1:N+1
    for j=1:1:M
    %Kinetic Energy
    f_n(i,1)=(1/2)*h(i,j).*(u(i,j).^2);
    lhs_integral = trapezoid(f_n,N,dx);
    Kinetic_Energy(i,j)=(rho*w)*lhs_integral;
    Kinetic_Energy(i,j)=Kinetic_Energy(i,j)/Kinetic_Energy(i,1);
    end
end
for i=1:1:N+1
    for j=1:1:M
    %Available Potential Energy
    f_n(i,1)=(1/2)*g*(eta(i,j).^2);
    lhs_integral = trapezoid(f_n,N,dx);
    Potential_Energy(i,j)=(rho*w)*lhs_integral;
    Potential_Energy(i,j)=Potential_Energy(i,j)/Potential_Energy(i,1);
    end
end
figure(5)
hold on
if(selection == 1)
    plot(t,Available_Energy)
    plot(t,Fluid_Mass)
    title('Available Energy and Fluid Mass')
end
if(selection == 0)
    plot(t,Available_Energy)
    plot(t,Fluid_Mass)
    title('Available Energy and Fluid Mass')

end
legend('Available Energy','Fluid Mass')
hold off

figure(6)
hold on
if selection ==1
    plot(t,Kinetic_Energy)
    plot(t,Available_Energy)
    plot(t,Potential_Energy)
    title('Kinetic, Available, and Potential Energy')
end

if selection ==0
    plot(t,Kinetic_Energy)
    plot(t,Available_Energy)
    plot(t,Potential_Energy)
    title('Kinetic, Available, and Potential Energy')
end
legend('Kinetic Energy', 'Available Energy', 'Potential Energy')
hold off
hold on
%run twice for linear and nonlinear subplots 
if selection==0
    figure(7)
    subplot(5,2,1);
    plot(x,h(:,1))
    title('m=1 Linear')
    xlabel('x');
    ylabel('h');
    subplot(5,2,3)
    plot(x,h(:,128))
    title('m=128 Linear')
    xlabel('x');
    ylabel('h');
    subplot(5,2,5)
    plot(x,h(:,256))
    title('m=256 Linear')
    xlabel('x');
    ylabel('h');
    subplot(5,2,7)
    plot(x,h(:,384))
    title('m=384 Linear')
    xlabel('x');
    ylabel('h');
    subplot(5,2,9)
    plot(x,h(:,512))
    title('m=384 Linear')
    xlabel('x');
    ylabel('h');
end
hold on   
if selection==1
    figure(7)
    subplot(5,2,2);
    plot(x,h(:,1))
    title('m=1 Non Linear')
    xlabel('x');
    ylabel('h');
    subplot(5,2,4)
    plot(x,h(:,128))
    title('m=128 Non Linear')
    xlabel('x');
    ylabel('h');
    subplot(5,2,6)
    plot(x,h(:,256))
    title('m=256 Non Linear')
    xlabel('x');
    ylabel('h');
    subplot(5,2,8)
    plot(x,h(:,384))
    title('m=384 Non Linear')
    xlabel('x');
    ylabel('h');
    subplot(5,2,10)
    plot(x,h(:,512))
    title('m=512 Non Linear')
    xlabel('x');
    ylabel('h');
end
hold off
