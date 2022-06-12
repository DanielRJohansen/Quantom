
%%
clear 
clc

data1 = readmatrix("D:\\Quantom\\potE.csv");
data2 = readmatrix("D:\\Quantom\\traj.csv");

%% RMSD

steps = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,24000,25000,26000,27000,28000,29000,30000,31000,32000,33000,34000,35000,36000,37000,38000,39000,40000,41000,42000,43000,44000,45000,46000,47000,48000,49000,50000,51000,52000,53000,54000,55000,56000,57000,58000,59000,60000,61000,62000,63000,64000,65000,66000,67000,68000,69000,70000,71000,72000,73000,74000,75000,76000,77000,78000,79000,80000,81000,82000,83000,84000,85000,86000,87000,88000,89000,90000,91000,92000,93000,94000,95000,96000,97000,98000,99000,100000,101000,102000,103000,104000,105000,106000,107000,108000,109000,110000,111000,112000,113000,114000,115000,116000,117000,118000,119000,120000,121000,122000,123000,124000,125000,126000,127000,128000,129000,130000,131000,132000,133000,134000,135000,136000,137000,138000,139000,140000,141000,142000,143000,144000,145000,146000,147000,148000,149000,150000,151000,152000,153000,154000,155000,156000,157000,158000,159000,160000,161000,162000,163000,164000,165000,166000,167000,168000,169000,170000,171000,172000,173000,174000,175000,176000,177000,178000,179000,180000,181000,182000,183000,184000,185000,186000,187000,188000,189000,190000,191000,192000,193000,194000,195000,196000,197000,198000,199000,200000,201000,202000,203000,204000,205000,206000,207000,208000,209000,210000,211000,212000,213000,214000,215000,216000,217000,218000,219000,220000,221000,222000,223000,224000,225000,226000,227000,228000,229000,230000,231000,232000,233000,234000,235000,236000,237000,238000,239000,240000,241000,242000,243000,244000,245000,246000,247000,248000,249000,250000,251000,252000,253000,254000,255000,256000,257000,258000,259000,260000,261000,262000,263000,264000,265000,266000,267000,268000,269000,270000,271000,272000,273000,274000,275000,276000,277000,278000,279000,280000,281000,282000,283000,284000,285000,286000,287000,288000,289000,290000,291000,292000,293000,294000,295000,296000,297000,298000,299000,300000,301000,302000,303000,304000,305000,306000,307000,308000,309000,310000,311000,312000,313000,314000,315000,316000,317000,318000,319000,320000,321000,322000,323000,324000,325000,326000,327000,328000,329000,330000,331000,332000,333000,334000,335000,336000,337000,338000,339000,340000,341000,342000,343000,344000,345000,346000,347000,348000,349000];
rmsd = [2.2821,6.81021,6.81333,6.74509,6.92436,6.74404,6.92011,6.91682,6.90913,6.90157,6.89176,6.92278,6.9113,6.86016,6.89209,6.88734,6.84267,6.83581,6.83203,6.83287,6.87386,6.87495,6.80278,6.81342,6.76878,6.72339,6.62601,6.65373,6.65036,6.64848,6.64829,6.64479,6.6435,6.64121,6.582,6.74603,6.81785,6.73266,6.72734,6.71907,6.64142,6.63576,6.67512,6.6625,6.71086,6.67992,6.73541,6.72469,6.70893,6.738,6.62788,6.67331,6.66137,6.54147,6.5236,6.55587,6.54684,6.5562,6.56823,6.656,6.47665,6.49826,6.29461,6.28139,6.11368,6.01665,5.93016,5.63033,5.36961,5.12214,4.99627,4.82671,4.68642,4.50096,4.34496,4.85623,5.53539,5.71856,6.05333,6.48002,7.42118,7.36134,8.12135,8.44244,8.72642,9.26082,9.03578,8.9387,8.76144,8.62922,8.57705,8.54548,8.42245,8.44438,8.32865,8.15725,8.19256,8.15326,8.14858,7.97877,7.64484,7.38017,7.11515,6.4951,6.58647,6.41539,6.6436,6.44001,6.26921,6.36219,6.48229,6.61705,6.68104,6.99016,7.01918,7.26949,7.25352,7.31652,7.16577,7.08625,7.30153,6.79211,7.16818,6.96476,6.84071,7.12816,7.64411,7.55033,7.42858,7.35461,7.52592,7.92201,7.77147,7.77318,7.69541,7.65802,7.73491,7.73587,7.79816,8.0078,8.16755,8.20765,8.29147,8.32131,8.34152,8.33457,8.21874,8.04551,8.098,7.99309,7.80217,7.59987,7.63652,7.55767,7.45286,7.37857,7.24104,6.97224,6.61915,6.44785,6.33047,6.20906,6.12919,6.0286,6.09893,6.1844,6.30736,6.41872,6.44871,6.58036,6.679,6.68227,6.72445,6.85195,6.70622,6.57592,6.58027,6.54089,6.69533,6.74481,6.60524,6.46067,6.42156,6.32038,6.11802,5.98455,5.86066,5.76732,5.66733,5.86338,5.72921,5.70408,6.0671,6.05436,6.39917,6.96918,7.15768,7.18239,7.34275,7.76138,7.80511,8.49198,8.5165,8.63802,8.63984,8.90306,9.02653,8.9913,8.93035,8.85156,8.52987,8.41964,8.52894,8.79725,8.655,8.72454,8.71407,8.445,8.25019,7.9077,8.06791,7.81778,7.62096,7.37465,7.21798,7.14659,7.55531,7.42686,6.87294,6.9654,6.95318,6.58703,6.0224,6.05646,5.79243,5.83025,5.42883,5.47846,4.89051,5.58154,6.51642,6.16412,5.27301,5.75043,6.42216,7.60799,7.55277,8.05927,7.9852,8.16355,7.67108,7.966,7.76516,7.92212,7.89716,7.5048,7.42003,7.67779,7.89835,7.94358,7.71462,7.56548,7.64396,7.86218,7.33801,7.34799,7.33541,7.29864,7.24925,7.19393,7.03787,6.92338,6.82747,6.75383,6.61508,6.4666,6.35922,6.16803,5.96186,5.74986,5.58444,5.41917,5.57826,5.4676,5.78399,5.9989,6.23221,5.63204,5.56055,6.14402,5.77154,5.6863,5.65203,5.69227,5.67465,5.68534,5.71926,5.37006,5.39312,4.64314,4.67265,4.2835,4.31738,4.3502,4.37927,4.26868,4.69656,4.66746,4.63043,4.77801,5.63722,5.57492,5.89233,5.82135,6.06585,6.77351,6.86613,6.96614,7.2168,7.63591,7.619,7.27718,7.04791,7.04155,7.14223,7.4638,6.81843,6.85266,6.16945,6.24086,5.6707,5.85751,5.52852,5.78194,5.77964,5.06471,5.24284,5.58186,5.59944,5.46882,5.46887,5.74082,5.99045,6.04049,6.1317,6.59127,7.0739,7.27269,7.39484];
plot(steps, rmsd, LineWidth=1.5)

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

workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_350000";
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
ylim([-inf 10000])

initial_temperature = kinE(1) * 2/(3*8.3145)
final_temperature = kinE(length(kinE)) * 2/(3*8.3145)

hold off


nexttile
t = (1:length(temp_data));
plot(t, temp_data)
%ylim([0 inf])
title("Temperature")
ylabel("Temperature [k]")
xlabel("time [fs]")



%% LOG data


clear 
clc

workdir = "C:\\PROJECTS\\Quantom\\Simulation\\Steps_350000";
file = fopen(strcat(workdir, "\\logdata.bin"), "rb");
data = fread(file, 'single');
fclose(file);

vel = data(1:10:end);
kin = data(2:10:end);
potE = data(3:10:end);
force = data(4:10:end);

x = 1:length(vel);
min_x = 340000;
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
%data = readtable("D:\\Quantom\\log.csv");
data = readtable("D:\\PROJECTS\\Quantom\\Simulation\\Steps_350000\\log.csv");

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

%% Torsion force
% CN1;NN2;CN9;HN9;0.000000;794.770020;3
phi = [-3.14:0.01:3.14];

%k = -7696.720215;
k = 794.77
ref = 0;

a = k*(1-cos(1*phi-ref));
b = k*(1-cos(2*phi-ref));
c = k*(1-cos(3*phi-ref));
%d = a+b+c;
%dd = -k*(sin(phi-ref) + 2.*sin(2.*phi-ref)*2 + 3.*sin(3.*phi-ref)*3);


%plot(phi, a)

%plot(phi, b)
%plot(phi, c)
plot(phi, a, LineWidth=1)
hold on
%plot(phi, a+b, LineWidth=1)
%plot(phi, a+b+c, LineWidth=1)

plot(phi, b, LineWidth=1)
plot(phi, c, LineWidth=1)

%plot(phi, dd)
%plot(phi, -k*sin(phi-ref))
ylabel("Potential Energy [J/mol]")
xlabel("Torsion [rad]")
legend("n=1", "n=2", "n=3")
%title("Torsion Potential")
grid 
hold off
%% LJ plot
clc
%sigma = 0.3923;
%epsilon = 0.5986*1000;

%Carbon:
sigma = 0.4
epsilon = 460
%x = 0.35:0.001:1;
x = 0.365:0.001:0.9;



y = 4 * epsilon * ((sigma./x).^12 - (sigma./x).^6);
dy = 24 .* epsilon .* 1./(x) .* (sigma./x).^6 .*(1-2.*(sigma./x).^6);

%tiledlayout(2,1)

%nexttile

plot(x,y, LineWidth=1.5)
title("Lennard-Jones potential between carbon atoms")
xlim([-inf, inf])
ylabel("Energy [J]")
xlabel("Distance [nm]")
grid

%% 

x = 0:0.01:1
sigma = 0.1
y = 1/sigma * exp(-0.5*((x-0.5)/sigma).^2)

px = [0.57 0.57]
py = [0 7.82]



plot(px, -py, LineWidth=2.5)
hold on
xlabel("Vector Activation")
ylabel("Loss")
plot(x,-y, LineWidth=1.5)
hold off




%%
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



