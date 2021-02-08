m_o = 9.109e-31;
m_n = 0.26*m_o;
T = 300;
kb = 1.380649e-23; %J/K
% 1 Electron Modelling
%1.1 vth = ?
v_th = sqrt(2*kb*T/m_n);

%1.2
tau_mn = 0.2e-12;
L_n = v_th*tau_mn;

%1.3

Xsize = 200e-9;
Ysize = 100e-9;
% Randomize position within region
p_x = rand(1000, 1)*Xsize;
p_y = rand(1000, 1)*Ysize;

% Give fixed velocity with random angle
v_x = zeros(1000, 1);
v_y = zeros(1000, 1);
theta = rand(1000, 1)*2*pi;
v_x = cos(theta).*v_th;
v_y = sin(theta).*v_th;

%Configure plot settings
line_color = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'];
clf
figure(1)
subplot(2, 1, 1)
xlim([0 Xsize])
ylim([0 Ysize])
title('Electron Movement through a Semiconductor')
xlabel('position (m)')
ylabel('position (m)')

delta_t = 5.35e-15; %time step size
time = linspace(0, 1000*delta_t, 1001);
T_avg = zeros(1001, 1);
%plot position then update
for i = 0:1000
    Xprev = p_x;
    Yprev = p_y;
    
    d_x = v_x .* delta_t;
    d_y = v_y .* delta_t;

    p_x = p_x + d_x;
    p_y = p_y + d_y;
    
    %Check boundary conditions
    for j = 1:1000
        %Horizontal boundaries
        if p_x(j) < 0
            p_x(j) = p_x(j) + Xsize;
        elseif p_x(j) > Xsize
            p_x(j) = p_x(j) - Xsize;
        end
        
        %Vertical boundaries
        if p_y(j) < 0
            v_y(j) = -v_y(j);
            p_y(j) = -p_y(j);
        elseif p_y(j) > Ysize
            v_y(j) = -v_y(j);
            p_y(j) = 2*Ysize-p_y(j);
        end
    end
    
    %plot points
    subplot(2, 1, 1)
    hold on
    for j = 1:14
        if abs(p_x(j)-Xprev(j)) <= v_th*delta_t
            %Only plot displacement if it is within maximum displacement,
            %otherwise it must have jumped across the horizontal boundary
            plot([p_x(j),Xprev(j)],[p_y(j),Yprev(j)], 'color', line_color(mod(j,7)+1))
        end
    end
    hold off
    
    %Measure temperature
    %E_k = 0.5*m*v^2 = 3/2*kb*T
    v = sqrt(v_x.^2 + v_y.^2);
    T_measured = (0.5.*m_n.*v.^2)./(kb.*3./2);
    T_avg(i+1) = mean(T_measured);
    subplot(2, 1, 2)
    plot(time(1:i+1), T_avg(1:i+1))
    title('Temperature of the Semiconductor')
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    pause(0.01)
end
hold off

%Part 2 will add the scattering of electrons

P_scat = 1 - exp(-delta_t/tau_mn);

% Randomize position within region
p_x = rand(1000, 1)*Xsize;
p_y = rand(1000, 1)*Ysize;

%Randomize velocities
theta = rand(1000, 1)*2*pi;
v_x = cos(theta).*v_th*randn();
v_y = sin(theta).*v_th*randn();

figure(2)
subplot(3, 1, 1)
xlim([0 Xsize])
ylim([0 Ysize])
title('Electron Scattering in a Semiconductor')
xlabel('position (m)')
ylabel('position (m)')

%plot position then update
for i = 0:1000
    Xprev = p_x;
    Yprev = p_y;
    
    for j = 1:1000
    	if P_scat > rand()
            %Particle scatters. Give new v_x and v_y
            v_x(j) = cos(2*pi*rand()).*v_th*randn();
            v_y(j) = sin(2*pi*rand()).*v_th*randn();
        end
    end
    
    d_x = v_x .* delta_t;
    d_y = v_y .* delta_t;
    p_x = p_x + d_x;
    p_y = p_y + d_y;
    
    %Check boundary conditions
    for j = 1:1000
        %Horizontal boundaries
        if p_x(j) < 0
            p_x(j) = p_x(j) + Xsize;
        elseif p_x(j) > Xsize
            p_x(j) = p_x(j) - Xsize;
        end
        
        %Vertical boundaries
        if p_y(j) < 0
            v_y(j) = -v_y(j);
            p_y(j) = -p_y(j);
        elseif p_y(j) > Ysize
            v_y(j) = -v_y(j);
            p_y(j) = 2*Ysize-p_y(j);
        end
    end
    
    v = sqrt(v_x.^2 + v_y.^2);
    subplot(3, 1, 3)
    hist(v, 20);
    title('Distribution of Electron Velocities')
    xlabel('Velocity (m/s)')
    ylabel('Number of Electrons')
    
    subplot(3, 1, 1)
    hold on 
    %plot points
    for j = 1:14
        if abs(p_x(j)-Xprev(j)) <= sqrt(2*(v_th^2))*delta_t
            %Only plot displacement if it is within maximum displacement,
            %otherwise it must have jumped across the horizontal boundary
            plot([p_x(j),Xprev(j)],[p_y(j),Yprev(j)], 'color', line_color(mod(j,7)+1))
        end
    end
    hold off
    
    
    %Measure temperature
    %E_k = 0.5*m*v^2 = 3/2*kb*T
    T_measured = (0.5.*m_n.*v.^2)./(kb.*3./2);
    T_avg(i+1) = mean(T_measured);
    subplot(3, 1, 2)
    plot(time(1:i+1), T_avg(1:i+1))
    title('Temperature of the Semiconductor')
    xlabel('Time (s)')
    ylabel('Temperature (K)')
    pause(0.01)
end


%Part 3 adds boxes of insulation that block the electron flow.

%box coordinates
box1_top = 100e-9;
box1_bottom = 60e-9;
box2_top = 40e-9;
box2_bottom = 0e-9;
box_left = 80e-9;
box_right = 120e-9;
top_box = [box_right, box1_bottom; box_right, box1_top; box_left, box1_top; box_left, box1_bottom; box_right, box1_bottom];
bottom_box = [box_right, box2_bottom; box_right, box2_top; box_left, box2_top; box_left, box2_bottom; box_right, box2_bottom];

% Randomize position within region
p_x = rand(1000, 1)*Xsize;
p_y = rand(1000, 1)*Ysize;

%Reposition any particles generated within the boxes
for i= 1:1000
    if (p_y(i) > box1_bottom) && (p_x(i) > box_left) && (p_x(i) < box_right)
        %particle is in top box
        p_x(i) = rand()*Xsize;
        p_y(i) = rand()*Ysize;
        %reduce i so that the updated particle position is re-tested
        i = i-1;
    elseif (p_y(i) < box2_top) && (p_x(i) > box_left) && (p_x(i) < box_right)
        %particle is in bottom box
        p_x(i) = rand()*Xsize;
        p_y(i) = rand()*Ysize;
        i = i-1;
    end
end

%Randomize velocities
theta = rand(1000, 1)*2*pi;
v_x = cos(theta).*v_th*randn();
v_y = sin(theta).*v_th*randn();

figure(3)
subplot(3, 1, 1)
xlim([0 Xsize])
ylim([0 Ysize])
title('Electron Movement through a Semiconductor with a Bottleneck')
xlabel('position (m)')
ylabel('position nm)')
hold on
%draw boxes
plot(top_box(:,1), top_box(:,2), 'color', 'k')
plot(bottom_box(:,1), bottom_box(:,2), 'color', 'k')
hold off

for i = 0:1000
    Xprev = p_x;
    Yprev = p_y;
    
    for j = 1:1000
    	if P_scat > rand()
            %Particle scatters. Give new v_x and v_y
            v_x(j) = cos(2*pi*rand()).*v_th*randn();
            v_y(j) = sin(2*pi*rand()).*v_th*randn();
        end
    end
    
    d_x = v_x .* delta_t;
    d_y = v_y .* delta_t;
    p_x = p_x + d_x;
    p_y = p_y + d_y;
    
    %Check boundary conditions
    for j = 1:1000
        %Horizontal boundaries
        if p_x(j) < 0
            p_x(j) = p_x(j) + Xsize;
        elseif p_x(j) > Xsize
            p_x(j) = p_x(j) - Xsize;
        end
        
        %Vertical boundaries
        if p_y(j) < 0
            v_y(j) = -v_y(j);
            p_y(j) = -p_y(j);
        elseif p_y(j) > Ysize
            v_y(j) = -v_y(j);
            p_y(j) = 2*Ysize-p_y(j);
        end
        
        %box boundaries
        if (p_y(j) > box1_bottom) && (p_x(j) > box_left) && (p_x(j) < box_right)
            %particle new position is in the top box
            if Xprev(j) > box_right && p_x(j) < box_right
                %right wall was hit
                v_x(j) = -v_x(j);
                p_x(j) = 2*box_right-p_x(j);
            elseif Xprev(j) < box_left && p_x(j) > box_left
                v_x(j) = -v_x(j);
                p_x(j) = 2*box_left-p_x(j);
            elseif Yprev(j) < box1_bottom && p_y(j) > box1_bottom
                v_y(j) = -v_y(j);
                p_y(j) = 2*box1_bottom-p_y(j);
            end
        elseif (p_y(j) < box2_top) && (p_x(j) > box_left) && (p_x(j) < box_right)
            %particle new position is in the bottom box
            if Xprev(j) > box_right && p_x(j) < box_right
                v_x(j) = -v_x(j);
                p_x(j) = 2*box_right-p_x(j);
            elseif Xprev(j) < box_left && p_x(j) > box_left
                v_x(j) = -v_x(j);
                p_x(j) = 2*box_left-p_x(j);
            elseif Yprev(j) > box2_top && p_y(j) < box2_top
                v_y(j) = -v_y(j);
                p_y(j) = 2*box2_top-p_y(j);
            end
        end
    end

    subplot(3, 1, 1)
    hold on 
    %plot points
    for j = 1:14
        if abs(p_x(j)-Xprev(j)) <= sqrt(2*(v_th^2))*delta_t
            %Only plot displacement if it is within maximum displacement,
            %otherwise it must have jumped across the horizontal boundary
            plot([p_x(j),Xprev(j)],[p_y(j),Yprev(j)], 'color', line_color(mod(j,7)+1))
        end
    end
    hold off
    
    pause(0.01)
end

density = zeros(100, 200);
v_sqr = zeros(100, 200);
%density and temperature
for i = 1:1000
    density(ceil(p_y(i)*1e9), ceil(p_x(i)*1e9)) = density(ceil(p_y(i)*1e9), ceil(p_x(i)*1e9)) + 1;
    v_sqr(ceil(p_y(i)*1e9), ceil(p_x(i)*1e9)) = v_sqr(ceil(p_y(i)*1e9), ceil(p_x(i)*1e9)) + v_x(i)^2 + v_y(i)^2;
end

subplot(3, 1, 2)
surf(density)
colormap(jet)
shading interp
title('Electron Density')
xlabel('position (nm)')
ylabel('position (nm)')
zlabel('Electrons per nm^2')

temp = m_n.*v_sqr./(2.*kb*3./2);
subplot(3, 1, 3)
surf(temp)
colormap(jet)
shading interp
title('Temperature')
xlabel('position (nm)')
ylabel('position (nm)')
zlabel('Temperature')
