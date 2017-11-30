clear all;

%1a

N      = [-1  
           1]; 

rmax = 1;  
Km = 100;

rates = @(rmax, Km, c)rmax * c(1) / (Km + c(1)) ; % c(1) = cS
       
dconc = @(rmax, Km, c) N * rates(rmax, Km, c) ;

odefun = @(t, c) dconc(rmax, Km, c);

tspan = 0:1000;    
c0 = [200 0];

[t_out,c_out] = ode45(odefun,tspan,c0');

figure('Name', 'Trajectories of states', 'Color', 'w');
plot(t_out, c_out);
legend({'S', 'P'});
xlabel('Time in sec');
ylabel ('Concentration in nM');

%1b

syms rmax cS Km cP

J_c = jacobian([((-rmax*cS)/(Km + cS)),((rmax*cS)/(Km + cS))],[cS, cP]) ;      % find jacobians in symbolic manner
J_p = jacobian([((-rmax*cS)/(Km + cS)),((rmax*cS)/(Km + cS))],[rmax, Km]) ;
     
clear rmax cS Km cP

rmax = 1;  
Km = 100;

%Jc = @(rmax, Km, c) J_c;
%Jp = @(rmax, Km, c) J_p;

Jc = @(rmax, Km, c) [-rmax/(Km+c(1))+rmax*c(1)/(Km+c(1))^2,0       % use of the above results in a non symbolic manner 
                     rmax/(Km+c(1))-rmax*c(1)/(Km+c(1))^2,0];     
                 
Jp = @(rmax, Km, c) [-c(1)/(Km+c(1)),+rmax*c(1)/(Km+c(1))^2
                      c(1)/(Km+c(1)),-rmax*c(1)/(Km+c(1))^2];

dsens = @(rmax,Km,c,s) Jc(rmax,Km,c)*s + Jp(rmax,Km,c);    

mat2vec = @(matrix) reshape(matrix',4,1);     
vec2mat = @(vector,rows) reshape(vector,rows,rows)';

odefun = @(t, y) [                                      % as given in the slides
    dconc(rmax, Km, y(1:2)) 
    mat2vec(dsens(rmax,Km,y(1:2),vec2mat(y(3:6),2))) 
];

s0 = [0 0 0 0];
[tout,c_s_out] = ode45(odefun,tspan,[c0 s0]');

figure('Name', 'Sensitivity', 'Color', 'w');
plot(tout, c_s_out(:,3:6));
legend({'dcS/drmax', 'dcP/drmax','dcS/dKm', 'dcP/dKm'});
xlabel('Time in sec');
ylabel ('Concentration in nM');

%1c

table1 = textread('dataset_1.txt');
hold on;
plot(table1(:,1),table1(:,2:3),'p') % plotted in points to be distinguished form the above plot
legend({'dcS/drmax', 'dcP/drmax','dcS/dKm', 'dcP/dKm','table1cS','table1cP'});

cov = [(2.5)^2 0
        0 25];
    
inv_cov = inv(cov);
F = zeros(2,2);

for i = 1:length(table1(:,1))
        S = vec2mat(c_s_out(table1(i, 1)+1,3:6),2); % time in each dataset corresponds to the index: time+1 in our c_s_out
        F = F + S' * inv_cov * S;
end

%1d

cr = sqrt(diag(inv(F)));
cr(1) = 100*cr(1)/rmax ;
cr(2) = 100*cr(2)/Km ;

disp(cr);

%1e

%%%%%%2nd table

table2 = textread('dataset_2.txt');
F2 = zeros(2,2);

for i = 1:length(table2(:,1))
        S = vec2mat(c_s_out(table2(i, 1)+1,3:6),2);
        F2 = F2 + S' * inv_cov * S;
end

cr2 = sqrt(diag(inv(F2)));
cr2(1) = 100*cr2(1)/rmax ; 
cr2(2) = 100*cr2(2)/Km ;

disp(cr2);

%%%%%%%%3rd table

table3 = textread('dataset_3.txt');

c0 = [100 0];
[tout, c_s_out] = ode45(odefun, tspan, [c0 s0]');

F3 = zeros(2,2);

for i = 1:length(table3(:,1))
        S = vec2mat(c_s_out(table3(i, 1)+1,3:6),2);
        F3 = F3 + S' * inv_cov * S;
end

cr3 = sqrt(diag(inv(F3)));
cr3(1) = 100*cr3(1)/rmax ;
cr3(2) = 100*cr3(2)/Km ;

disp(cr3);
