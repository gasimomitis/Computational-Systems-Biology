clear all;

%1a

N      = [-1  1  0  0  0  0  0  0  0  0  % stoichiometric matrix
           1 -1  0  0  0  0  0  0  0  0 
           0  0 -1  0  0  1  0  0  0  0 
           0  0  1 -1  1 -1  0  0  0  0 
           0  0  0  1 -1  0  0  0  0  0 
           0  0  0  0  0  0 -1  0  0  1 
           0  0  0  0  0  0  1 -1  1 -1
           0  0  0  0  0  0  0  1 -1  0]; 

V = [2.5;0.25;0;0;0.75;0.75;0;0;0.5;0.5];  % values as denoted in tables
k = [0;0;0.025;0.025;0;0;0.025;0.025;0;0];
K = [10;8;15;15;15;15;15;15;15;15];

rates = @(V,K,k,c)[V(1)*c(1)/(K(1)+c(1))   % like in the toy model
                   V(2)*c(2)/(K(2)+c(2))
                   k(3)*c(2)*c(3)/(K(3)+c(3))
                   k(4)*c(2)*c(4)/(K(4)+c(4))
                   V(5)*c(5)/(K(5)+c(5))
                   V(6)*c(4)/(K(6)+c(4))
                   k(7)*c(5)*c(6)/(K(7)+c(6))
                   k(8)*c(5)*c(7)/(K(8)+c(7))
                   V(9)*c(8)/(K(9)+c(8))
                   V(10)*c(7)/(K(10)+c(7))] ;
       
dconc = @(V, K, k, c) N * rates(V,K,k,c) ;

%1b

odefun = @(t, c) dconc(V, K, k, c);

tspan = 1:1000;     % concetrations are computed in 1000 time values in sec
c0 = [100 0 300 0 0 300 0 0]';

[t_out,c_out] = ode45(odefun,tspan,c0);

figure('Name', 'Trajectories of all states', 'Color', 'w');
plot(t_out, c_out);
legend({'MKKK', 'MKKK_P', 'MKK', 'MKK_P', 'MKK_{PP}', 'MAPK', 'MAPK_P', 'MAPK_{PP}'});
xlabel('Time in sec');
ylabel ('Concentration in nM');

%1c

threshold = 0.000001;  % threshold value difference between concecutive c_out
f_time = 100; % time length
tspan = 1:f_time;
V1 = 0:0.02:0.5 ; 
results = zeros(length(V1),1);

for i = 1:length(V1)
	V(1) = V1(i); 
    %dconc = @(V, K, k, c) N * rates(V,K,k,c) ;
    odefun = @(t, c) dconc(V, K, k, c);
    c0 = [100 0 300 0 0 300 0 0]';
    flag = true;
    while (flag)
        [t_out,c_out] = ode45(odefun,tspan,c0);
        flag = ~(abs(c_out(f_time,:)-c_out(f_time-1,:)) < threshold) ;  % check if concentration at the last two time moments differ
        c0 = c_out(f_time,:)' ; % set initial condition for next time group the current concetrations. 
    end
    results(i) = c_out(f_time,8);
end

figure('Name', 'Steady state values of MAPK_{PP} as a function of the V1', 'Color', 'w');
plot(V1, results);
xlabel('V1');
ylabel ('MAPK_{PP} Concentration in Steady Stade');


