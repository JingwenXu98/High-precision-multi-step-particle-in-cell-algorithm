%Boris
clear
clc
close all
format long

q = 1;
m = 1; 
Grid = 1e8;
out = 10000;
dt = 1e-3;
T = Grid*dt;
t = 0:dt:T;

x = zeros(3 , Grid+1);
v = zeros(3 , Grid+1);

E=[0; 0; 0]; 
B=[0; 0; -0.2];

x(:,1)=[0; 0; 0] ; 
v(:,1)=[2; 0; 0];
tic

%boris coefficients
Bx=[0,B(3),-B(2);-B(3),0,B(1);B(2),-B(1),0];
Ek(1,:)=0.5 * sum((v(1,:).^2));

for i = 1:Grid
%step1
   v1 = v(:,i) + (q*E*dt)/(2*m); %v1=v^{-}
   v2 = (eye(3) - q*dt/(2*m) * Bx) \ (eye(3) + q*dt/(2*m) * Bx) * v1 ; %v2=v^{+}
   v(:,i+1) = v2 + (q*E*dt)/(2*m);
%step2  
   v(:,i+1) = v(:,i) + dt * q/m * (E + Bx * (v(:,i+1)+v(:,i))/2);
   x(:,i+1) = x(:,i) + dt * v(:,i+1);
%Ene   
   Ek(i+1,:)=0.5 * sum(v(:,i+1).^2);
   if mod(i,out) ==0
       i
       toc
   end
end

%% figure
figure(1)
plot3(x(1,:),x(2,:),x(3,:),'m')
hold on

figure(2)
plot3(v(1,:),v(2,:),v(3,:),'m')
hold on

figure(3)
plot(Ek,'m')
hold on