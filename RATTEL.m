%     RATTEL-RAT Timidity Tester Evidence based Learning-To find optimal
%     location of deployment of TOMCAT-our deployment device.
%     Copyright (C) 2018  Yash Rana,Devang Liya,Niteshwar M.A.
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

clear;

a=150;%Dimensions of the arena.a is the width, b is the heigth.
b=100;
m=1;%Mass of the rats.

%kb=;%Boltzman's Constant.
%T=;% Effective Temperature-Can model the acitivity of the rat.
%R=;%Radius-Will model the size of the rat.
%n=;%Viscosity of the medium.
%D=((k*T)/(6*n*R*pi));%Diffusion Constant;

Dt=2;% Translational DIffusion Constant-For Now!
Dr=0.001;%Rotational Diffusion Constant-For Now!

%Using some approximate values of the parameters.
%We make the box (or the arena);


delx=10;%Step size for Fearomone Field
nx=a/delx;%No. of steps in x.
ny=b/delx;%No. of steps in y.

ts=0.0;%Time (s) for start of supply of Fearomone
tc=100;%Time step for cutoff of supply of Fearomone

R=0.5;%Rate of production of Fearomone in TOMCAT.

C=zeros(nx,ny,2);%Fearomone Field.

Dp=20;%For now.Diffusion rate of Fearomone.

%T= input('Enter the Time Period: ');     % Time (in Seconds)
T=10;
N=100.*T;% Number of Impulses/change in track (it's just fancy word for 'steps')
h=sqrt(T/N);%T/N is the time step;
%np=input('Enter the Number Of Particle: ');    % Number of Particles
np=500;

nf=input('Enter the Number Of Grain Stacks: ');    % Number of Grain Stacks

fl=zeros(nf,2); %Location of food grain stacks.

for i=1:nf
  fl(i,:)=input(['Enter the coordinates of Grain Stack No. ' num2str(i) ' in the form [x,y] (0<x<150) (0<y<100) : ']);    
end
    
% Initialization of the position of particles.

x = zeros(1,np); %Let's go with this logic.
y = zeros(1,np);
theta=zeros(1,np);

%We need to define the initial point for the different rats and also their
%individual orientations.
 
%Here we make a function to decide the initial position.This function
%starts of all the rats along the walls with their orientation facing away
%from the wall

[x,y]=INI(x,y,np,a,b,theta);

av(1,:)=1*ones(1,np);%Active velocity component
vx=zeros(1,np);%X component of interaction velocity.
vy=zeros(1,np);%X component of interaction velocity.

% Iteration to store positions of particles

rng('shuffle'); %To re-initialise the random number generator.

%The following code segment is to compare all the various possible locations where TOMCAT can be deployed.

for m=2:nx-1 
  for n=2:ny-1

C=zeros(nx,ny,2);

for i=1:N
    o(1,i)=0;%Order Parmeter.Will be calculated for each location of the TOMCAT.
  
%For each time step we first make the pheromone field.

 for Fi=2:nx-1
     for Fj=2:ny-1
         
         if((Fi==m) & (Fj==n) & (i>=(ts.*100)) & (i<=(tc.*100)))
      
             C(Fi,Fj,2)=(Dp/(delx*delx))*(C(Fi+1,Fj,1)+C(Fi-1,Fj,1)+C(Fi,Fj+1,1)+C(Fi,Fj-1,1)-4*C(Fi,Fj,1))+C(Fi,Fj,1)+R;
            
         else 
             C(Fi,Fj,2)=(Dp/(delx*delx))*(C(Fi+1,Fj,1)+C(Fi-1,Fj,1)+C(Fi,Fj+1,1)+C(Fi,Fj-1,1)-4*C(Fi,Fj,1))+C(Fi,Fj,1);
            
         end
         
     end
 end
 
 %Now we work on the boundaries.
 
%These are Nueman Boundary Conditions-Partial Reflective Boundary Condition. 
  C(1,:,2)=C(2,:,2)-(0.005*delx);%We work around the boundaries in an anticlockwise fashion.
  C(:,1,2)=C(:,2,2)-(0.005*delx);
  C(nx,:,2)=C(nx-1,:,2)-(0.005*delx);
  C(:,ny,2)=C(:,ny-1,2)-(0.005*delx);
 
  %Now we update the coordinates of all the active rats.
for j=1:np

x(i+1,j)=x(i,j)+av(1,j)*cos(theta(i,j))*h*h+h*sqrt(2*Dt)*randn()+vx(i,j)*h*h; 


%The randn function generates arrays of random numbers whose elements are normally distributed with mean 0, variance 1, and standard deviation 1. 
y(i+1,j)=y(i,j)+(av(1,j)*sin(theta(i,j))*h*h)+(h*sqrt(2*Dt)*randn())+vy(i,j)*h*h;



%Now we add the reflection conditions.
                                                    
if y(i+1,j)<0
    y(i+1,j)=-1*y(i+1,j);

else if (y(i+1,j)>(b))
    y(i+1,j)=((2*b)-y(i+1,j));
else if (x(i+1,j)<0)
     x(i+1,j)=-1*x(i+1,j);
else if (x(i+1,j)>(a))
     x(i+1,j)=((2*a)-x(i+1,j));               
    end
    end
    end
end

%We repeat the reflections condition to take care of the rare case of
%multiple reflections.

                              
if y(i+1,j)<0
    y(i+1,j)=-1*y(i+1,j);

else if (y(i+1,j)>(b))
    y(i+1,j)=((2*b)-y(i+1,j));
else if (x(i+1,j)<0)
     x(i+1,j)=-1*x(i+1,j);
else if (x(i+1,j)>(a))
     x(i+1,j)=((2*a)-x(i+1,j));               
    end
    end
    end
end


%We need to determine the force acting on the particles.The

%Fx(i,j),Fy(i,j),depend on x(i,j) and y(i,j).These change v(i+1);
%We will include all the forces: fear,attraction to food,repulsion from the
%pheromone in this Force calculation.

[Fx(i,j),Fy(i,j)]=Force(x(i,j),y(i,j),C,delx,nx,ny,nf,fl); 

%%SINCE WE ARE FOLLOWING THE ACTIVE-BROWNIAN PARTICLE MODEL,WE CAN'T UPDATE
%%VELOCITIES VIA THIS METHOD-True!

%Now we update the "Interaction" velocities.

vx(i+1,j)=vx(i,j)+((Fx(i,j)/m)*h*h);
vy(i+1,j)=vy(i,j)+((Fy(i,j)/m)*h*h);

%Now we update the theta.No tourque applied as of now.

theta(i+1,j)=theta(i,j)+h*sqrt(2*Dr)*randn();%Will include memory effects in future.

%Now we calculate the Order Parameter.

for u=1:nf
    o(1,i)=o(1,i)+sqrt(  (x(i,j)-fl(u,1))^2  +(y(i,j)-fl(u,2))^2 );
end
    
end

o(1,i)=o(1,i)/(nf*np*180.27);% where (nf*np*180.27) is the normalization constant.

%Resetting the Pheromone field.Very essential.
C(:,:,1)=C(:,:,2);

end

FO(m-1,n-1)=mean(o);

  end
end

figure;
surf(FO');%The was we imagine space is the transpose of the way it is organized in memroy.
view(2);
axis([1 13 1 8]);
colorbar;
figure;
maxMatrix = max(FO(:));
[row,col] = find(FO==maxMatrix);
row=row+1;
col=col+1;
disp(['Hence the most optimum location is [' num2str(row*10) ',' num2str(col*10) ']'] );


%Now we would like to make the video only of the run where the order parameter
%is the maximum.
 

V = VideoWriter('20181017_Rat_Model_FINAL_CHECK_1'); % Video object.
V.FrameRate=2;
open(V); 

[Y,X] = meshgrid(0:delx:b-1,0:delx:a-1);
 
hfig=figure('Visible','off');
fid=figure;

%For the purpose of preseting the dynamics to the user we make np=50;
np=50;
cmap = hsv(np); % Creates a np-by-3 set of colors from the HSV colormap
x = zeros(1,np); %Let's go with this logic.
y = zeros(1,np);
theta=zeros(1,np);
%We need to define the initial point for the different rats and also their
%individual orientations.
 %Here we make a function to decide the initial position.
[x,y,theta]=INI(x,y,np,a,b,theta);

av=1*ones(1,np);%Active velocity component
vx=zeros(1,np);
vy=zeros(1,np);
C=zeros(nx,ny,2);

for i=1:N


 %For each time step we first make the pheromone field. 
  for Fi=2:nx-1
      for Fj=2:ny-1
          
          if((Fi==row) & (Fj==col) & (i>=(ts.*100)) & (i<=(tc.*100)))
       
              C(Fi,Fj,2)=(Dp/(delx*delx))*(C(Fi+1,Fj,1)+C(Fi-1,Fj,1)+C(Fi,Fj+1,1)+C(Fi,Fj-1,1)-4*C(Fi,Fj,1))+C(Fi,Fj,1)+R;
             
          else 
              C(Fi,Fj,2)=(Dp/(delx*delx))*(C(Fi+1,Fj,1)+C(Fi-1,Fj,1)+C(Fi,Fj+1,1)+C(Fi,Fj-1,1)-4*C(Fi,Fj,1))+C(Fi,Fj,1);
             
          end
          
      end
  end
  
%Now we work on the boundaries.
 % These are Nueman Boundary Conditions-Partial reflective Boundary Condition. 
  C(1,:,2)=C(2,:,2)-(0.005*delx);%We work around the boundaries in an anticlockwise fashion.
 C(:,1,2)=C(:,2,2)-(0.005*delx);
 C(nx,:,2)=C(nx-1,:,2)-(0.005*delx);
 C(:,ny,2)=C(:,ny-1,2)-(0.005*delx);
 
 
for j=1:np

x(i+1,j)=x(i,j)+av(1,j)*cos(theta(i,j))*h*h+h*sqrt(2*Dt)*randn()+vx(i,j)*h*h; 


%The randn function generates arrays of random numbers whose elements are normally distributed with mean 0, variance 1, and standard deviation 1. 
y(i+1,j)=y(i,j)+(av(1,j)*sin(theta(i,j))*h*h)+(h*sqrt(2*Dt)*randn())+vy(i,j)*h*h;



%Now we add the reflection conditions.
                                                    
if y(i+1,j)<0
    y(i+1,j)=-1*y(i+1,j);

else if (y(i+1,j)>(b))
    y(i+1,j)=((2*b)-y(i+1,j));
else if (x(i+1,j)<0)
     x(i+1,j)=-1*x(i+1,j);
else if (x(i+1,j)>(a))
     x(i+1,j)=((2*a)-x(i+1,j));               
    end
    end
    end
end

%We repeat the reflections condition to take care of the rare case of
%multiple reflections.

                              
if y(i+1,j)<0
    y(i+1,j)=-1*y(i+1,j);

else if (y(i+1,j)>(b))
    y(i+1,j)=((2*b)-y(i+1,j));
else if (x(i+1,j)<0)
     x(i+1,j)=-1*x(i+1,j);
else if (x(i+1,j)>(a))
     x(i+1,j)=((2*a)-x(i+1,j));               
    end
    end
    end
end


%We need to determine the force acting on the particles.The

%Fx(i,j),Fy(i,j),depend on x(i,j) and y(i,j).These change v(i+1);
%We will include all the forces: fear,attraction to food,repulsion from the
%pheromone in this Force calculation.
[Fx(i,j),Fy(i,j)]=Force(x(i,j),y(i,j),C,delx,nx,ny,nf,fl); 

%%SINCE WE ARE FOLLOWING THE ACTIVE-BROWNIAN PARTICLE MODEL,WE CAN'T UPDATE
%%VELOCITIES VIA THIS METHOD-True!

%Now we update the "Interaction" velocities.

vx(i+1,j)=vx(i,j)+((Fx(i,j)/m)*h*h);
vy(i+1,j)=vy(i,j)+((Fy(i,j)/m)*h*h);

%Now we update the theta.No tourque applied as of now.

theta(i+1,j)=theta(i,j)+h*sqrt(2*Dr)*randn();%Will include memory effects in future.

    
end
%Resetting the Pheromone field.Very essential.
C(:,:,1)=C(:,:,2);

if(rem(i,10)==0)

    ax1 = axes('Position',[0.1 0.1 0.7 0.7]);%The big graph for the rats
    ax2 = axes('Position',[0.64 0.64 0.28 0.28]);%The graph to show the diffusion.
    fid.CurrentAxes=ax1;
    axis equal;
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');

%Now we have the coordinates of the rats for the best location of the TOM-CAT. 

%figure('Visible','off');
    rectangle(ax1,'Position',[0 0 a b],'EdgeColor','k','LineWidth',0.1); 
%axes(ax1);
    hold on;
    grid on;
%axes(ax1);
    axis(ax1,[-5 a+5 -5 b+5]);
    text(ax1,10,10,num2str(i));
    scatter(ax1,fl(:,1),fl(:,2),100,'green','filled');  
%Drawing the TOM-CAT Pheromone Box. 
    scatter(ax1,row*10,col*10,150,'red','filled');
 
    for k=1:np
    scatter(ax1,x(i,k),y(i,k),15,cmap(k,:),'filled');
    %Drawing tail 
    if(i<21)
        plot(ax1,x(1:i,k),y(1:i,k),'Color',cmap(k,:));
    else
        plot(ax1,x(i-20:i,k),y(i-20:i,k),'Color',cmap(k,:));
    end

    end
  
  
    %For visualization of the field. We would like to see the Fearomone field.
    figure(fid);
    fid.CurrentAxes=ax2;
    surf(ax2,X,Y,C(:,:,2),'FaceAlpha',0.5);
    axis(ax2,[0 150 0 100 0 5]);
    
    text(ax2,5,5,0.1,num2str(k));
 
    frame=getframe(gcf);
    size(frame.cdata);
    writeVideo(V,frame);
    clf;
 end
  


end

close(V);

function Po=Pot(x,y,nf,fl) %Function to define the potential in which the rats move.

a=150;%Dimensions of the arena.a is the width, b is the heigth.
b=100;

A=750;%Parameter for the Gaussian defining thpotential for the fear of open spaces and love for food.
B=800;

%Let's unleash the KRAKEN!

Po=A*exp((-((x-(a/2))^2)-((y-(b/2))^2))/B);

for i=1:nf
  Po=Po-1*A*exp((-((x-(fl(i,1)))^2)-((y-(fl(i,2)))^2))/B);
end

end

function [Fx,Fy]=Force(x,y,C,delx,nx,ny,nf,fl) %Function to calculate force from Potential acting on the active rats.

dx=0.01;%dx for calculating gradient of potential.
dy=0.01;

c1=5;%Scaling the Static force.

c2=5000;%Scaling of the dynamic forces.


%Now we need to calculate the gradient of the potential to find the force.

tsx=-1*((Pot(x+dx,y,nf,fl)-Pot(x-dx,y,nf,fl))/(2*dx)); %All these are for static potentials; 
tsy=-1*((Pot(x,y+dy,nf,fl)-Pot(x,y-dy,nf,fl))/(2*dy));

%Now we include the dynamic potential of the pheromone.

%Let us first find the box in the nx*ny boxes in which our rat lies.This
%will also be useful when we introduce memory.

Ri=floor(x/delx) +1;
Rj=floor(y/delx) +1;

if (Ri==1) 
    
    tdx=-1*((C(Ri+1,Rj,2)-C(Ri,Rj,2))/delx);

else if (Ri==nx)
        
    tdx=-1*((C(Ri,Rj,2)-C(Ri-1,Rj,2))/delx);
    
    else
        
        tdx=-1*((C(Ri+1,Rj,2)-C(Ri-1,Rj,2))/(2*delx));
    end
end

if (Rj==1)
    
    tdy=-1*((C(Ri,Rj+1,2)-C(Ri,Rj,2))/delx);

else if (Rj==ny)
        
    tdy=-1*((C(Ri,Rj,2)-C(Ri,Rj-1,2))/delx);
    
    else
        
        tdy=-1*((C(Ri,Rj+1,2)-C(Ri,Rj-1,2))/(2*delx));
    end
end

Fx=c1*tsx+c2*tdy;
Fy=c1*tsy+c2*tdy;

end

function [x,y,theta]=INI(x,y,np,a,b,theta)%This function is calculate the initial locations of the rats along the walls along with the their orientations.
 c=randi([0 3],1,np);

 for i=1:np
     if(c(i)==0)
       x(1,i)=a*rand();
       y(1,i)=b-1;
       
       %vx(1,i)=rand(-5,5);
       %vy(1,i)=rand(-5,0);
       theta(1,i)=pi+(pi*rand());
     
     else if (c(i)==1)
             x(1,i)=1;
             y(1,i)=b*rand();
             %vx(1,i)=rand(0,5);
             %vy(1,i)=rand(-5,5);
             theta(1,i)=-1*(pi/2)+pi*rand();
         else if (c(i)==2)
             x(1,i)=a*rand();
             y(1,i)=1;
             %vx(1,i)=rand(-5,5);
             %vy(1,i)=rand(0,5);
             theta(1,i)=pi*rand();
             
             else if(c(i)==3)
                 x(1,i)=a-1;
                 y(1,i)=b*rand();
                 %vx(1,i)=rand(-5,0);
                 %vy(1,i)=rand(-5,5);
                  theta(1,i)=(pi/2)+(pi*rand());
             
                 end
             end
         end
     end
 end
 
end
