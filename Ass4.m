
clear
clc
%% Question 1 + Question 2 
% I do not want mistakes from my incomplete Assingmnet 3 to affect this to
% carry over. R3 was set to 10. 

 
%% Question 3

R1 = 1;
R2 = 2;
C = 0.25;
L = 0.2;
alpha = 100;
R4 = 0.1;
R3 = 10; 
Ro = 1000;

G = [
    1   0   0   0   0   0   0;
    -1/R2   (1/R1)+(1/R2)   -1   0   0   0   0;
    0   1   0   -1   0   0   0;
    0   0   -1   1/R3   0   0   0;
    0   0   0   0   -alpha   1   0;
    0   0   0   1/R3   -1   0   0;
    0   0   0   0   0   -1/R4   (1/R4)+(1/Ro)];

C_Mat = [
    0   0   0   0   0   0   0;
    -C   C   0   0   0   0   0;
    0   0   -L   0   0   0   0;
    0   0   0   0   0   0   0;
    0   0   0   0   0   0   0;
    0   0   0   0   0   0   0;
    0   0   0   0   0   0   0];

% V = [V2; V1; Il; V3; I3; V4; Vo]

%i DC Sweep 
Vin = linspace(-10,10,100);
Vo = zeros(100,1);
V3 = zeros(100,1); 

for i = 1:100
    F =[Vin(i); 0; 0; 0; 0; 0; 0];
    V = G\F;
    Vo(i) = V(7);
    V3(i) = V(4);
end

figure(1)
plot(Vin,Vo)
title('Vin vs Vout')
xlabel('Vin (V)')
ylabel('Vout(V)')

figure(2)
plot(Vin,V3)
title('Vin vs V3')
xlabel('Vin (V)')
ylabel('V3 (V)')


%ii AC Sweep

omega = linspace(1,1E4,10000);
Vin = 1; 
Vo = zeros(10000,1);
 F =[Vin; 0; 0; 0; 0; 0; 0];
for i = 1:1000
    V = (G +(omega(i)*2*pi*C_Mat *1j)) \F;
    Vo(i) = V(7);
end

figure(3)
semilogx(omega,Vo)
title('Omega vs Vout')
xlabel('omega')
ylabel('Vout(V)')

gain = 20*log10(abs(Vo/Vin));
figure(4)
semilogx(omega,gain)
title('Omega vs gain')
xlabel('Vin (V)')
ylabel('Vout(V)')

%iii AC Sweep2 (with random perturbations on C)

omeega = pi;

for i = 1:10000
    C_dist = normrnd(C,0.05);
    C_Mat = [
        0   0   0   0   0   0   0;
        -C_dist   C_dist   0   0   0   0   0;
        0   0   -L   0   0   0   0;
        0   0   0   0   0   0   0;
        0   0   0   0   0   0   0;
        0   0   0   0   0   0   0;
        0   0   0   0   0   0   0];
    
    V = (G +(omega(i)*2*pi*C_Mat *1j)) \F;
    Vo(i) = V(7);
end

gain = 20*log(abs(Vo/Vin));
figure(5)
histogram(gain)
title('Gain with Capacitor Noise')
xlabel('Gain (dB)');

%% Question 4a
% An RLC Circuit

%% Question 4b
% The capacitance and iductance in the cicuit are affected by the
% frequency. As the frequency is increased, the impedance of the cpacitor
% (Zc) decreases, while the impedance of the Inductor (Zl) Increases.

%% Question 4d
% Finite difference solution in the Time Domain: 

simulation_time = 1000; % run simulation for 1 second 
dt = 1; % time step 

inputA = zeros(1,1000);
Vpre = zeros(7,1); %V(i-1)
% start simulation 
for time = dt:dt:simulation_time
   
    % input A: step that transitions from 0 to 1 at t=0.03s
    if time < 30
        input_A(1,time) = 0;
    else 
        input_A(1,time) = 1; 
    end
    
     F =[input_A(1,time); 0; 0; 0; 0; 0; 0];
     
     Vout = (G +C)\ (C.*Vpre) + F;
     Vout1_store(1,time) = Vout(7);
     
     Vpre = Vout;
end

Time = linspace(dt,simulation_time,simulation_time);

figure(6)
plot(Time,input_A)
hold on;
plot(Time,abs(Vout1_store))
xlim([0 1000])
ylim([0 1.2])
xlabel('Time (ms)')
ylabel('Voltage(V)')
hold off
legend('Vin','Vout')

  % input B: sin(2?f t) function with f = 1/(0.03) 
input_B = zeros(1,1000);
Vpre = zeros(7,1); %V(i-1)
f = 1/30;
   
% start simulation 
for time = dt:dt:simulation_time
   
     input_B(time) = sin(2*pi*f*time);
    
     F =[input_B(1,time); 0; 0; 0; 0; 0; 0];
     
     Vout = (G +C)\ (C.*Vpre) + F;
     Vout2_store(1,time) = Vout(7);
     
     Vpre = Vout;
end

Time = linspace(dt,simulation_time,simulation_time);

figure(7)
plot(Time,input_B)
hold on;
plot(Time,abs(Vout2_store))
xlim([0 100])
ylim([-1.1 1.1])
xlabel('Time (ms)')
ylabel('Voltage(V)')
hold off
legend('Vin','Vout')


    
% input C: A guassian pulse with a magnitude of 1, std dev. of 0.03s
% and a delay of 0.06s
input_C = zeros(1,1000);
   
% start simulation 
for time = dt:dt:simulation_time
   
    input_C(1,time) = (exp(-2*log(2)*(time-60).^2/(30)^2)).*cos(-2*pi*f*(time-60));
    
     F =[input_C(1,time); 0; 0; 0; 0; 0; 0];
     
     Vout = (G +C)\ (C.*Vpre) + F;
     Vout3_store(1,time) = Vout(7);
     
     Vpre = Vout;
end

Time = linspace(dt,simulation_time,simulation_time);

figure(8)
plot(Time,input_C)
hold on;
 plot(Time,abs(Vout3_store))
  xlim([0 200])
%  ylim([-1.1 1.1])
xlabel('Time (ms)')
ylabel('Voltage(V)')
hold off
legend('Vin','Vout')


% Frequncy Plots 
A_fin = fftshift(fft(input_A));
B_fin = fftshift(fft(input_B)); 
C_fin = fftshift(fft(input_C));

A_fout = fftshift(fft(Vout1_store));
B_fout = fftshift(fft(Vout2_store));
C_fout = fftshift(fft(Vout3_store));

figure(9)
plot(Time,A_fin)
hold on
plot(Time,A_fout)
hold off
title('Fourier Transform: Step Input')
xlabel('Frequency')
ylabel('Magnitude')
legend('Input','Output')

figure(10)
plot(Time,B_fin)
hold on
plot(Time,B_fout)
hold off
title('Fourier Transform: Sine Wave Input')
xlabel('Frequency')
ylabel('Magnitude')
legend('Input','Output')

figure(11)
plot(Time,C_fin)
hold on
plot(Time,C_fout)
hold off
title('Fourier Transform: Gaussian Pulse Input')
xlabel('Frequency')
ylabel('Magnitude')
legend('Input','Output')


%% Question 5

Cn = 0.00001
In = 0,001 % initial I0
G = [
    1   0   0   0   0   0   0;
    -1/R2   (1/R1)+(1/R2)   -1   0   0   0   0;
    0   1   0   -1   0   0   0;
    0   0   -1   1/R3   In   0   0;
    0   0   0   0   -alpha   1   0;
    0   0   0   1/R3   -1   0   0;
    0   0   0   0   0   -1/R4   (1/R4)+(1/Ro)];

C_Mat = [
    0   0   0   0   0   0   0;
    -C   C   0   0   0   0   0;
    0   0   -L   0   0   0   0;
    0   0   0   Cn   0   0   0;
    0   0   0   0   0   0   0;
    0   0   0   0   0   0   0;
    0   0   0   0   0   0   0];

% V = [V2; V1; Il; V3; I3; V4; Vo]

for time = dt:dt:simulation_time
   
    input_C(1,time) = (exp(-2*log(2)*(time-60).^2/(30)^2)).*cos(-2*pi*f*(time-60));
    
     F =[input_C(1,time); 0; 0; 0; 0; 0; 0];
     
     Vout = (G +C)\ (C.*Vpre) + F;
     Vout3_store(1,time) = Vout(7);
     
     Vpre = Vout;
     In = normrnd(0.001,0.01);
end

Time = linspace(dt,simulation_time,simulation_time);

figure(12)
plot(Time,input_C)
hold on;
plot(Time,abs(Vout3_store))
  xlim([0 200])
%  ylim([-1.1 1.1])
title('Vout and Vin w/ Cn and random In')
xlabel('Time (ms)')
ylabel('Voltage(V)')
hold off
legend('Vin','Vout')

figure(13)
plot(Time,fft(Vout3_store))
%  ylim([-1.1 1.1])
title('Fourier Transform of Output')
xlabel('Time (ms)')
ylabel('Voltage(V)')
hold off
dt = 100;
for time = 1:dt:simulation_time
   
    input_C(1,time) = (exp(-2*log(2)*(time-60).^2/(30)^2)).*cos(-2*pi*f*(time-60));
    
     F =[input_C(1,time); 0; 0; 0; 0; 0; 0];
     
     Vout = (G +C)\ (C.*Vpre) + F;
     Vout3_store(1,time) = Vout(7);
     
     Vpre = Vout;
     In = normrnd(0.001,0.01);
end

Time = linspace(dt,simulation_time,simulation_time);

figure(14)
plot(Time,input_C)
hold on;
plot(Time,abs(Vout3_store))
  xlim([0 200])
%  ylim([-1.1 1.1])
title('Vout and Vin w/ Cn and random In (1/100 time step')
xlabel('Time (ms)')
ylabel('Voltage(V)')
hold off
legend('Vin','Vout')

figure(15)
plot(Time,fft(Vout3_store))
%  ylim([-1.1 1.1])
title('Fourier Transform of Output (1/100 time step)')
xlabel('Time (ms)')
ylabel('Voltage(V)')
hold off
dt = 10;
for time = 1:dt:simulation_time
   
    input_C(1,time) = (exp(-2*log(2)*(time-60).^2/(30)^2)).*cos(-2*pi*f*(time-60));
    
     F =[input_C(1,time); 0; 0; 0; 0; 0; 0];
     
     Vout = (G +C)\ (C.*Vpre) + F;
     Vout3_store(1,time) = Vout(7);
     
     Vpre = Vout;
     In = normrnd(0.001,0.01);
end

Time = linspace(dt,simulation_time,simulation_time);

figure(16)
plot(Time,input_C)
hold on;
plot(Time,abs(Vout3_store))
  xlim([0 200])
%  ylim([-1.1 1.1])
title('Vout and Vin w/ Cn and random In (1/10s step')
xlabel('Time (ms)')
ylabel('Voltage(V)')
hold off
legend('Vin','Vout')

figure(17)
plot(Time,fft(Vout3_store))
%  ylim([-1.1 1.1])
title('Fourier Transform of Output (1/10s Time step)')
xlabel('Time (ms)')
ylabel('Voltage(V)')
hold off

