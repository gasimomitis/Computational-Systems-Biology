clear all

% a)

%Two state model: K1, K1P

N = [-1 0  1 0     % stoichiometric matrix
      1 0 -1 0];
     
    %K1;K2;Feed;P;P_mass
k = [30;30;0;3;0.75*10^8];

    %K1;K2;Feed;P
K = [10^(-6);10^(-6);0;2.5*10^(-9)];

rates = @(K,k,c)[ k(1) * c(3) .* c(1) ./ (Km(1) + c(1))
                  k(2) * c(2) .* c(1) ./ (Km(2) + c(1))
                  k(4) * c(4) .* c(3) ./ (Km(4) + c(3))
                  k(5) * c(4) .* c(3)] ;
              
dconc = @(K, k, c) N * rates(K,k,c) ; 

% b)

c_K1_T = 5*10^(-8);
c_P = 5*10^(-9);
c_K1P = 0:10^(-9):c_K1_T;
c_K1 = c_K1_T - c_K1P;

r_K1 = (k(1) * c_K1P .* c_K1) ./ (K(1) + c_K1) ;
r_P = (k(4) * c_P .* c_K1P) ./ (K(4) + c_K1P) ;

figure('Name', 'r_K1, r_P', 'Color', 'w');
plot(c_K1P,r_K1,'Color','b');
hold on
plot(c_K1P,r_P,'Color','r');
xlabel('concetration K_{1P}');
legend({'r_{K1}', 'r_{P}'})

% c)

clear c_K1P c_K1 r_K1 r_P

syms c_K1P c_K1
dc_K1 = sym('-(k_K1 * c_K1P * c_K1) / (K_K1 + c_K1) + (k_P * c_P * c_K1P) / (K_P + c_K1P)');
dc_K1P = sym ('k_K1 * c_K1P * c_K1 / (K_K1 + c_K1) - (k_P * c_P * c_K1P) / (K_P + c_K1P)');
nullcline_K1 = solve(dc_K1 ,c_K1);    % we have solved for C_K1 in both nullclines
nullcline_K1P = solve(dc_K1P ,c_K1);

 % we have solved for C_K1 in both nullclines. Nullclines appear to be same. In
 % order to find the ss we need to find the points in which the two
 % nullclines intersect. So we need one more equation. This is the
 % conservation of the total concetration of K1. C_K1 = C_K1_T - C_K1P.

ss_K1P   = solve(sym('c_K1_T - c_K1P') - nullcline_K1, c_K1P);

syms ss_K1

ss_K1_1  = solve(sym('c_K1_T - ss_K1') - ss_K1P(1),ss_K1);
ss_K1_2  = solve(sym('c_K1_T - ss_K1') - ss_K1P(2),ss_K1);


ss_K1_eval1 = eval(subs(ss_K1_1,{'k_K1','k_P','K_K1','K_P','c_P','c_K1_T'}, [k(1) k(4) K(1) K(4) c_P c_K1_T]));
ss_K1_eval2 = eval(subs(ss_K1_2,{'k_K1','k_P','K_K1','K_P','c_P','c_K1_T'}, [k(1) k(4) K(1) K(4) c_P c_K1_T]));
ss_K1P_eval1 = eval(subs(ss_K1P(1),{'k_K1','k_P','K_K1','K_P','c_P','c_K1_T'}, [k(1) k(4) K(1) K(4) c_P c_K1_T]));
ss_K1P_eval2 = eval(subs(ss_K1P(2),{'k_K1','k_P','K_K1','K_P','c_P','c_K1_T'}, [k(1) k(4) K(1) K(4) c_P c_K1_T]));

figure('Name', 'Steady States', 'Color', 'w');
plot([ss_K1P_eval1 ss_K1P_eval2],[ss_K1_eval1 ss_K1_eval2],'p');
xlabel('concetration K_{1P}');
ylabel('concetration K_{1}');

% d)

syms c_K1P c_K1
dc_K1 = sym('-(k_K1 * c_K1P * c_K1) / (K_K1 + c_K1) + (k_P2 * c_P * c_K1P)');
dc_K1P = sym ('k_K1 * c_K1P * c_K1 / (K_K1 + c_K1) - (k_P2 * c_P * c_K1P)');
nullcline_K1 = solve(dc_K1 ,c_K1);    % we have solved for C_K1 in both nullclines
nullcline_K1P = solve(dc_K1P ,c_K1);

 % we have solved for C_K1 in both nullclines. Nullclines are the same. In
 % order to find the ss we need to find the points in which the two
 % nullclines intersect. So we need one more equation. This is the
 % conservation of the total concetration of K1. C_K1 = C_K1_T - C_K1P.

ss_K1P   = solve(sym('c_K1_T - c_K1P') - nullcline_K1, c_K1P);

syms ss_K1

ss_K1_1  = solve(sym('c_K1_T - ss_K1') - ss_K1P,ss_K1);


ss_K1_eval = eval(subs(ss_K1_1,{'k_K1','k_P2','K_K1','c_P','c_K1_T'}, [k(1) k(5) K(1) c_P c_K1_T]));
ss_K1P_eval = eval(subs(ss_K1P(1),{'k_K1','k_P2','K_K1','c_P','c_K1_T'}, [k(1) k(5) K(1) c_P c_K1_T]));

figure('Name', 'Steady States', 'Color', 'w');
plot(ss_K1P_eval,ss_K1_eval,'p');
xlabel('concetration K_{1P}');
ylabel('concetration K_{1}');

% we observe 1 ss less



