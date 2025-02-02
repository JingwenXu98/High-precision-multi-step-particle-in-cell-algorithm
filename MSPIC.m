clear;clc;close all;format long;

dt=1e-3;
t=0:dt:1e+3;
N=length(t)-1;

x=zeros(N+1,3);
v=zeros(N+1,3);
E=zeros(N+1,3);
B=zeros(N+1,3);


w_E=[0,0,0];
w_B=[0,0,-0.2];


B(:,3) =  w_B(3);
v(1,1) = 2;
 
q=1;
m=1;

%Class RK
for n=1:3
%     %%%%%sin()
%     E_RK_1 = sin( w_E .*(n+1/2)* dt) ;
%     E_RK_2 = sin( w_E .*(n+1)* dt) ;
%     B_RK_1 = sin( w_B .*(n+1/2)* dt) ;
%     B_RK_2 = sin( w_B .*(n+1)* dt) ;
%     %%%%%cos()
%     E_RK_1 = 2 * cos( w_E .*(n+1/2)* dt) ;
%     E_RK_1(2:3)=0;
%     E_RK_2 = 2 * cos( w_E .*(n+1)* dt) ;
%     E_RK_2(2:3)=0;
%     B_RK_1 = 2 * cos( w_B .*(n+1/2)* dt) ;
%     B_RK_2 = 2 * cos( w_B .*(n+1)* dt) ;

%     %%%%Constant
    E_RK_1 = w_E ;
    E_RK_2 = w_E ;
    
    B_RK_1 = w_B ;
    B_RK_2 = w_B ;

%     %%%%%Linear
%     E_RK_1 = w_E .*(n+1/2)* dt;
%     E_RK_2 = w_E .*(n+1)* dt ;
%     
%     B_RK_1 = w_B ;
%     B_RK_2 = w_B ;
    
    Q_1 = v(n,:);
    k_1 = q/m * ( E(n,:)+ cross(Q_1,B(n,:)) );
    
    Q_2 = v(n,:) + k_1 * dt/2;
    k_2 = q/m * ( E_RK_1+ cross(Q_2,B_RK_1) );
    
    Q_3 = v(n,:) + k_2 * dt/2;
    k_3 = q/m * ( E_RK_1+ cross(Q_3,B_RK_1) );
    
    Q_4 = v(n,:) + k_3 * dt;
    k_4 = q/m * ( E_RK_2+ cross(Q_4,B_RK_2) );
  
    
    v(n+1,:) = v(n,:) + dt/6 * ( k_1 + 2*k_2 + 2*k_3 + k_4 );
    x(n+1,:) = x(n,:) + dt/6 * ( Q_1 + 2*Q_2 + 2*Q_3 + Q_4 );
    
end
%%
tic
for n= 4:N       
%%%%%%%%Adams-Bashforth-Hamming%%%%%%%%%%%%%%%%%%%  
v(n+1,:) = v(n,:) + dt/24 * ( -9 * q / m * ( E(n-3,:) + cross(v(n-3,:),B(n-3,:)) ) + 37 * q / m * ( E(n-2,:) + cross(v(n-2,:),B(n-2,:)) ) ...
                              - 59 * q / m * ( E(n-1,:) + cross(v(n-1,:),B(n-1,:)) )  + 55 * q / m * ( E(n,:) + cross(v(n,:),B(n,:)) )  );
                          
x(n+1,:) = 9/8 * x(n,:) - 1/8 * x(n-2,:) + dt * 3/8 * (- v(n-1,:) + 2 * v(n ,:) + v(n+1,:) );           
                          
v(n+1,:) = 9/8 * v(n,:) - 1/8 * v(n-2,:) + dt * 3/8 * ( - q / m * ( E(n-1,:) + cross(v(n-1,:),B(n-1,:)) ) + 2 * q / m * ( E(n,:) + cross(v(n,:),B(n,:)) ) ...
                              + q / m * ( E(n+1,:) + cross(v(n+1,:),B(n+1,:)) ) );  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

      if (mod(n,5000)==0) %show time every 5000 steps
        toc
        disp(['time: ',num2str(n)])
      end                         
end  

Ene= 1/2 * (v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
%% Section : Plot figure
figure(1)
plot3(v(:,1),v(:,2),v(:,3))
hold on
figure(2)
plot3(x(:,1),x(:,2),x(:,3))
hold on
figure(3)
plot(Ene)