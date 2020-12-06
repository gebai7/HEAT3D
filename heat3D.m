clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The simulation of three-dimensional heat conduction equations is used 
% to learn and understand such problems, including two types of boundary 
% conditions that can be used: neuman and Dirichlet boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%input
density = 3000;  % kg/m^3
conductivity = 350;   % W/(mK)
heat_capacity = 900;  % J/(kgK)


alpha = conductivity/(density*heat_capacity);

NN=1000;
L_x = 10;
L_y = 10;
L_z = 10;
N_x = 20;
N_y = 20;
N_z = 20;
dx = L_x/N_x;
dy = L_y/N_y;
dz = L_z/N_z;
dt_ideal = 1/(2*alpha*(1/dx^2+1/dy^2+dz^2));
% dt_ideal1 = 1/(2*alpha1*(1/dx^2+1/dy^2));
% dt_ideal2 = 1/(2*alpha2*(1/dx^2+1/dy^2));
fprintf('Critical time-step size : %.2f\n\n',dt_ideal);

dt = 500;  % CHANGE TIME-STEP SIZE HERE %
% Initial condition and boundary conditions
%Dirichlet boundary
%%%notice : If you want to use Neuman boundary conditions
%%% do not need to annotate the following boundary
T_initial = 100;
T_top =100;    % 
T_bottom = 300; % 
T_left =150;    % 
T_before = 200; % 
T_behind =100;    % 
T_right = 50;  % 


% Applying initial and boundary conditions to the grid
T = zeros(N_x+1,N_y+1,N_z+1,NN);
T(:,:,1,:) = T_bottom;
T(:,:,N_z+1,:) = T_top;
T(1,:,:,:) = T_left;
T(N_x+1,:,:,:) = T_right;
T(:,1,:,:) = T_behind;
T(:,N_y+1,:,:) =T_before;
T(:,:,:,1)=T_initial;
% alpha_z=zeros(N_x+2,N_y+2);
% alpha_z(:,:)=alpha;
% alpha_z(40:end,:)=alpha2;
% alpha_z(40:60,35:45)=alpha1;
%% Steady Solution %%
% Euler Method
k = 1;
    error_max(k) = 100;
    error_min(k) = 100;
for   k=1:NN
        for i = 2:N_x
            for j = 2:N_y
                for q = 2:N_z
                    T(i,j,q,k+1) = T(i,j,q,k)+dt*alpha*(((T(i-1,j,q,k)-2*T(i,j,q,k)+T(i+1,j,q,k))/dx^2)+((T(i,j-1,q,k)-2*T(i,j,q,k)+T(i,j+1,q,k))/dy^2)+((T(i,j,q-1,k)-2*T(i,j,q,k)+T(i,j,q+1,k))/dz^2));
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%Neuman boundary
%         for i = 2:N_x
%             for j = 2:N_y
%                 for q = 2:N_z
%                     T(i,1,q,k+1)= (T(i,2,q,k+1)+T(i,3,q,k+1))/2;%ºó
%                     T(1,j,q,k+1)= (T(2,j,q,k+1)+T(3,j,q,k+1))/2;%×ó
%                     T(i,N_y+1,q,k+1)= (T(i,N_y,q,k+1)+T(i,N_y-1,q,k+1))/2;%Ç°
%                     T(N_x+1,j,q,k+1)= (T(N_x,j,q,k+1)+T(N_x-1,j,q,k+1))/2;%ÓÒ
%                 end
%             end
%         end

end
% T=T(:,:,:,1:k);

sol_time = k*dt;
% T_min = 100;
% T_max = 600;
T_max = max([T_top,T_bottom,T_right,T_left,T_initial,T_before,T_behind]);
T_min = min([T_top,T_bottom,T_right,T_left,T_initial,T_before,T_behind]);

%% ------------------------------ SOLUTION ----------------------------- %%


%-------------------------------------------------------------------------%

%% Generating the plate %%
[x,y,z]=meshgrid(1:N_x+1,1:N_y+1,1:N_z+1);
tt=T(:,:,:,1);
sx =15:20;sy=15:20;sz=1:5;

%% -------------------------- Plotting Section ------------------------- %%
for l=1:10:NN
figure(5)
% subplot(2,2,3)
% slice(x,y,z,T(:,:,:,l))
slice(x,y,z,T(:,:,:,l),sx,sy,sz);colormap('jet');
alpha()
grid off
shading interp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end