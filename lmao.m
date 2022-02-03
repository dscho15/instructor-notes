clear; clc;
n = 3;

syms q1 q2 q3 qd1 qd2 qd3 a1 a2 a3 alpha1 alpha2 alpha3 d1 d2 d3 real

% DH matrix
DH = [a1, alpha1, d1, q1, qd1; ...
      a2, alpha2, d2, q2, qd2; ...
      a3, alpha3, d3, q3, qd3];

% Inertia matrices
I1 = zeros(3);
I1(1,1) = 0.0084;
I1(2,2) = 0.0064;
I1(3,3) = 0.0084;

I2 = zeros(3);
I2(1,1) = 0.0078;
I2(2,2) = 0.21;
I2(3,3) = 0.21;

I3 = zeros(3);
I3(1,1) = 0.0016;
I3(2,2) = 0.0462;
I3(3,3) = 0.0462;

I = cell(n, 1);
I{1} = I1;
I{2} = I2;
I{3} = I3;

% DH parameters
%DH = zeros(3,5);
%DH(1,2) = pi/2;
%DH(1,3) = 0.1625 
%DH(2,1) = -0.425;
%DH(3,1) = -0.3922;
%DH(:, 4) = q;
%DH(:, 5) = qd;

% Mass
m = [3.761,8.058,2.846];

p_com = cell(n, 1);
p_com{1} = [0,-0.02561,0.00193]';
p_com{2} = [0.2125,0,0.11336]';
p_com{3} = [0.15,0,0.0265]';

g = [0, 0, -9.82]';

%%

L = Lagrange(DH, m, p_com, I, g);

L_a = subs(L, [a1 a2 a3 alpha1 alpha2 alpha3 d1 d2 d3], ...
    [0, -0.425, -0.3922, pi/2, 0, 0, 0.1625, 0, 0]);

dLdq = gradient(L_a, [q1 q2 q3]);
dLdqd = gradient(L_a, [qd1 qd2 qd3]);

syms q1_f(t) q2_f(t) q3_f(t) 
qd1_f(t) = diff(q1_f, t);
qd2_f(t) = diff(q2_f, t);
qd3_f(t) = diff(q3_f, t);

dLdq_f = subs(dLdq, [q1, q2, q3, qd1, qd2, qd3], [q1_f, q2_f, q3_f, qd1_f, qd2_f, qd3_f]);
dLdqd_f = diff(subs(dLdqd, [q1, q2, q3, qd1, qd2, qd3], [q1_f, q2_f, q3_f, qd1_f, qd2_f, qd3_f]), t);

Q = dLdqd_f - dLdq_f;

syms qdd1 qdd2 qdd3 real

%f = subs(Q, [diff(q1_f(t), t, t), diff(q2_f(t), t, t), diff(q3_f(t), t, t)], [qdd1, qdd2, qdd3]);
%f = subs(f, [diff(q1_f(t), t), diff(q2_f(t), t), diff(q3_f(t), t)], [qd1, qd2, qd3]);
%f = subs(f, [q1_f(t), q2_f(t), q3_f(t)], [q1, q2, q3]);
%f = vpa(f, 4);

%%

%vpa(subs(f, [q1 q2 q3 qd1 qd2 qd3 qdd1 qdd2 qdd3], [1 pi/3 pi/3 0 0 0 0 0 0]), 4)

%%

D = -eye(3);


qdd1_v = solve(Q == [0; 0; 0], [diff(q1_f(t), t, t); diff(q2_f(t), t, t); diff(q3_f(t), t, t)])


%%

tspan = [0 5];
y0 = 0;
[t,y] = ode45(f, tspan, [1 pi/3 pi/3 0 0 0 0 0 0]);

%%

function b_p_com = CoM_of_link(b_T_link, link_p_com)
    b_R_link = b_T_link(1:3, 1:3);
    b_p_com = b_T_link(1:3, 4) + b_R_link * link_p_com;
end

function [Jp, Jo, R] = jac_to_CoM(DH, link_p_com, i)
    
    % Transformation to all link frames
    A = sym(eye(4));
    Ai = cell(i+1, 1);
    Ai{1} = A;
    for idx = 1:i
        A = A * DH2trans(DH(idx, 1), DH(idx, 2), DH(idx, 3), DH(idx, 4));
        Ai{idx+1} = A;
    end
    
    b_p_com = CoM_of_link(Ai{i+1}, link_p_com);
    
    % Initialize Jacobian
    Jp = sym(zeros(3, size(DH, 1)));
    Jo = sym(zeros(3, size(DH, 1)));
    
    for idx = 1:i
        A = Ai{idx};
        z = A(1:3, 3);
        p = A(1:3, 4);
        
        Jp(:, idx) = cross(z, b_p_com - p);
        Jo(:, idx) = z;
    end
    
    R = Ai{i+1}(1:3, 1:3);

end

function A = DH2trans(a, alpha, d, theta)

    % Translation along Z
    Tz = sym(diag(ones(1, 4)));
    Tz(3, 4) = d;
    
    % Rotation about Z
    Rz = [cos(theta), -sin(theta), 0, 0;...
          sin(theta), cos(theta), 0, 0; ...
          0, 0, 1, 0;...
          0, 0, 0, 1];
    
    % Translation along X
    Tx = sym(diag(ones(1, 4)));
    Tx(1, 4) = a;
    
     % Rotation about X
    Rx = [1, 0, 0, 0; ...
          0, cos(alpha), -sin(alpha), 0; ...
          0, sin(alpha), cos(alpha), 0; ...
          0, 0, 0, 1];
    
    A = Rz * Tz * Rx * Tx;
end

function E = potential(DH, g, m, p_com)
    
    n = size(DH, 1);
    
    E = 0;
    A = eye(4);
    for i = 1:n
        
        A = A * DH2trans(DH(i, 1), DH(i, 2), DH(i, 3), DH(i, 4));
        p = CoM_of_link(A, p_com{i});
        
        E = E - m(i) * g' * p;
    end
end

function E = kinetic(DH, m, p_com, I)
    n = size(DH, 1);
    
    qdot = DH(:, 5);
    
    E = 0;
    A = eye(4);
    for i = 1:n
        
        [Jp, Jo, R] = jac_to_CoM(DH, p_com{i}, i);
        
        A = A * DH2trans(DH(i, 1), DH(i, 2), DH(i, 3), DH(i, 4));
        p = CoM_of_link(A, p_com{i});
        
        Il = R * I{i} * R';
        
        E = E + 0.5 * m(i) * qdot' * Jp' * Jp * qdot ...
            + 0.5 * qdot' * Jo' * Il * Jo *  qdot; 
    end
end

function L = Lagrange(DH, m, p_com, I, g)

    L = kinetic(DH, m, p_com, I) - potential(DH, g, m, p_com);
    
end
