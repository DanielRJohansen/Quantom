
%%
clear 
clc

data1 = readmatrix("D:\\Quantom\\potE.csv");
data2 = readmatrix("D:\\Quantom\\traj.csv");

%%
key_particle = 2;

potE = data1(:,key_particle);
traj = data2(:,key_particle*3-2:key_particle*3);
x = 1:length(potE);

v = traj(2:length(x),:) - traj(1:length(x)-1, :);
v = v.^2;
v = sum(v,2);
v = sqrt(v);
v(v>1) = 0;
x_v = x(2:length(x));




tiledlayout(2,1)
nexttile
plot(x, potE)
title('Total LJ Potential (repulsion) on particle')
ylim([-4000 200000])

nexttile
plot(x_v,v)
ylim([0 0.01])

%% Force changes in 1 dt
x=categorical([1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2097152 4194304 8388608 16777216 33554432 67108864 134217728 268435456 536870912 ]);
y=[48842 921 1438 2246 2669 1552 547 391 352 303 170 94 82 55 39 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];

bar(x,log10(y))

%% Energy revision
clear 
clc

workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_6000";
file = fopen(strcat(workdir, "\\energy.bin"), "rb");
energy_data = fread(file, 'single');
fclose(file);

file = fopen(strcat(workdir, "\\temperature.bin"));
temp_data = fread(file, 'single');
fclose(file);

length(energy_data)
n_elements = length(energy_data)/3
energy_data = reshape(energy_data, [3, n_elements])';
%size(data);

%t = (1:length(temp_data)).*200;
%plot(t, temp_data)
%ylim([-2000 2000])
 


tiledlayout(2,1)
nexttile
potE = energy_data(:,1);
kinE = energy_data(:,2);
totalE = energy_data(:,3);
%totalE = kinE + potE;
x = 1:length(potE);



plot(x, potE);
hold on
plot(x, kinE);
plot(x, totalE);
title("Average energy")
legend("Potential energy", "Kinetic energy", "Total energy");
ylabel("Energy [J/mol]")
xlabel("time [fs]")
%xlim([0 10000])
%label("Kinetic energy")

%xlim([00 2700])
%ylim([000 20000])

initial_temperature = kinE(1) * 2/(3*8.3145)
final_temperature = kinE(length(kinE)) * 2/(3*8.3145)

hold off


nexttile
t = (1:length(temp_data));
plot(t, temp_data)
ylim([0 inf])
title("Temperature")
ylabel("Temperature [k]")
xlabel("time [fs]")



%% LOG data


clear 
clc

workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_200000";
file = fopen(strcat(workdir, "\\logdata.bin"), "rb");
data = fread(file, 'single');
fclose(file);

vel = data(1:10:end);
kin = data(2:10:end);
potE = data(3:10:end);
force = data(4:10:end);

x = 1:length(vel);
min_x = 000;
max_x = length(x);


tiledlayout(3,1)
nexttile
plot(x, vel)
title('Velocity on particle')
ylabel("[nm/ns]")
xlabel("[fs]")
xlim([min_x max_x])

nexttile
plot(x, potE)
hold on
plot(x, kin)
plot(x, kin+potE)
title('Energy in particle')
ylabel("[J/mol]")
xlabel("[fs]")
xlim([min_x max_x])
legend("Potential energy", "Kinetic energy", "Total energy");

hold off




nexttile
plot(x, force)
title('Force')
ylabel("[N]")
xlabel("[fs]")
xlim([min_x max_x])

%% Waterforce
clear 
clc
data = readtable("D:\\Quantom\\log.csv");

%tiledlayout(3,1)
%nexttile
f_x = data.Data8(2:end);
f_y = data.Data9(2:end);
f_z = data.Data10(2:end);

t = 1:length(f_x);
plot(t, f_x);
hold on
plot(t, f_y);
plot(t, f_z);
title("Waterforce on single particle")
legend("Force x", "Force y", "Force z");
ylabel("Force [N/mol]")
xlabel("time [fs]")
hold off

%% LJ plot
clc
sigma = 0.3923;
epsilon = 0.5986*1000;
%x = 0.35:0.001:1;
x = 0.244:0.001:1;



y = 4 * epsilon * ((sigma./x).^12 - (sigma./x).^6);
dy = 24 .* epsilon .* 1./(x) .* (sigma./x).^6 .*(1-2.*(sigma./x).^6);

tiledlayout(2,1)

nexttile
plot(x,y)
ylabel("[J]")
xlabel("[nm]")

nexttile
plot(x,dy)
ylabel("[N]")
xlabel("[nm]")


%% Test vibrations
clear
clc




len = 2000;
t = 0:len;
dt = 30 * 10^-6;

acc = zeros(1,2);
vel = zeros(1,2);
force = zeros(1,2);
force_prev = zeros(1,2);

mass = [15.999 1.008 ] ./ 1000;
pos = [0 0.099];

pos_history = zeros(length(t), 2);
bond_len = zeros(length(t));




force(1) = calcForce(pos(1), pos(2));
force(2) = calcForce(pos(2), pos(1));
for i = 1:len+1
    pos(1) = integratePos(pos(1), dt, vel(1), mass(1), force(1));
    pos(2) = integratePos(pos(2), dt, vel(2), mass(2), force(2));
    force_prev = force;


    force(1) = calcForce(pos(1), pos(2));
    force(2) = calcForce(pos(2), pos(1));
    pos;
    vel(1) = integrateVel(vel(1), dt, mass(1), force(1), force_prev(1));
    vel(2) = integrateVel(vel(2), dt, mass(2), force(2), force_prev(2));


    pos_history(i, 1) = pos(1);
    pos_history(i, 2) = pos(2);
    bond_len(i) = pos(2)-pos(1);

end

tiledlayout(2,1)
nexttile
plot(t, bond_len);
title("Bond len, 1D")
ylabel("Bond len [nm]")
xlabel("Step #")
%ylim([0.09 inf])


nexttile
plot(t, pos_history(:,2));
title("atom positions");
hold on
plot(t, pos_history(:,1));



function [pos] = integratePos(pos_, dt, vel, mass, force)
pos = pos_ + dt * (vel + 0.5/mass * force * dt);
end

function [vel] = integrateVel(vel_, dt, mass, force, force_prev)
vel = vel_ + dt * 0.5 / mass * (force + force_prev);
end





function [acc, vel, pos] = integrate(acc_prev, vel_prev, pos_prev, force, mass)
    dt = 2.5 * 10^-6;
    pos = pos_prev + vel_prev * dt + acc_prev * 0.5 * dt * dt;

    acc_next = force * (1000./mass);
    vel_next = vel_prev + (acc_prev + acc_next) * 0.5 * dt;

    acc = acc_next;
    vel = vel_next;

end

function [y] = norm(x)
y = x/abs(x);
end

function [force] = calcForce(pos0, pos1) %repulsive
kb = 17.5 * 10e+6;
ref = 0.095;
dist = abs(pos0-pos1);
dif = dist - ref;

invert = 1;
if dif > 0
    invert = -1;
end

scalar= 0.5 * kb * (dif*dif);
force = norm(pos0-pos1) * scalar * invert;
end


