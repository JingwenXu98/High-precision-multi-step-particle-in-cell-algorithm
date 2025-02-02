%DirectIm
clear;clc;close all;format long;

Deltat = [1e-1, 9e-2, 8e-2, 7e-2, 6e-2, 5e-2, 4e-2, 3e-2, 2e-2,...
          1e-2, 9e-3, 8e-3, 7e-3, 6e-3, 5e-3, 4e-3, 3e-3, 2e-3,... 
          1e-3, 9e-4, 8e-4, 7e-4, 6e-4, 5e-4, 4e-4, 3e-4, 2e-4,...
          1e-4, 9e-5, 8e-5, 7e-5, 6e-5, 5e-5, 4e-5, 3e-5, 2e-5,... 
          1e-5, 9e-6, 8e-6, 7e-6, 6e-6, 5e-6, 4e-6, 3e-6, 2e-6,...
          1e-6, 9e-7, 8e-7, 7e-7, 6e-7, 5e-7, 4e-7, 3e-7, 2e-7,...
          1e-7];
EneDirectIm = zeros(1,55);
q = 1;
m = 1; 
Grid = 1e+8;
out = 1e6;

for n = 19:19
    
dt = Deltat(n);
T = Grid*dt;
t = 0:dt:T;

x = zeros(3 , Grid+1);
v = zeros(3 , Grid+1);
A = zeros(3 , Grid+1);

E=[0; 0; 0]; 
B=[0; 0; -0.2];

x(:,1)=[0; 0; 0]; 
v(:,1)=[2; 0; 0];

tic
EDIm(1,:)=0.5 * sum((v(1,:).^2));

%<T>
Omega = q * B * dt /2 /m;
TB = [1 + Omega(1)^2, Omega(1) * Omega(2) + Omega(3),  Omega(1) * Omega(3) - Omega(2);...
      Omega(1) * Omega(2) - Omega(3), 1 + Omega(2)^2,  Omega(2) * Omega(3) + Omega(1);...
      Omega(1) * Omega(3) + Omega(2), Omega(2) * Omega(3) - Omega(1), 1 + Omega(3)^2 ] ./ (1 + sum(Omega.^2)) ;
A(:,1) = q * E /m;
for i = 1:Grid
    %Step 1
    A(:,i+1) = (A(:,i) + q * E /m)/2;
    v(:,i+1) = TB * (v(:,i) + A(:,i) * dt /2 + cross(v(:,i),Omega));
    x(:,i+1) = x(:,i) + v(:,i+1) * dt;
    %Step 2
    v(:,i+1) = v(:,i+1) + q/2/m * TB * E *dt;
    x(:,i+1) = x(:,i+1) + q/2/m * E *dt^2;
    
    EDIm(i+1,:)=0.5 * sum(v(:,i+1).^2);
   if mod(i,out) ==0
       i
       toc
   end
end

%% Ene
currentFile = sprintf('EneDirectIm%d.mat',dt);
save(currentFile,'EDIm')
EneDirectIm(n) = max(abs(EDIm-EDIm(1)));

end
currentFile = sprintf('AllEneDt.mat');
save(currentFile,'EneDirectIm')
%%
figure(1)
plot3(x(1,:),x(2,:),x(3,:))

figure(2)
plot3(v(1,:),v(2,:),v(3,:))

figure(3)
plot(EDIm)