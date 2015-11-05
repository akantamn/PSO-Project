function [L_Js, x_bar, Lambda_bar] = calculatelagrangian(a, b, c, w, u, P,Pdiff,timestep, iter)
% INPUTS: u(size m), a(m), b(m), c(m), w(T,m+1), T
% OUTPUTS: L_Js % L_Js(j,:) = calculatelagrangian(a, b, c, w, u, T);


% Pmin = [50, 100, 100];
% Pmax = [600, 400, 200];
m = size(a,2); % number of variables
n = size(w,1); % number of constraints


%     L(j) = (a(1)*x(j,1)*x(j,1) + b(1)*x(j,1) + c(1))*u(j,1) ...
%     + (a(2)*x(j,2)*x(j,2) + b(2)*x(j,2) + c(2))*u(j,2) ...
%     + (a(3)*x(j,3)*x(j,3) + b(3)*x(j,3) + c(3))*u(j,3) ...
%     + Lambda(j)*(w(j,4) + w(j,1)*x(j,1)*u(j,1) + w(j,2)*x(j,2)*u(j,2)+ w(j,3)*x(j,3)*u(j,3));

%     dLdx1 = (2*a(1)*x(1) + b(1) + Lambda*w(j,1)) * u(1) = 0
%     dLdx2 = (2*a(2)*x(2) + b(2) + Lambda*w(j,2)) * u(2) = 0
%     dLdx3 = (2*a(3)*x(3) + b(3) + Lambda*w(j,3)) * u(3) = 0
%     dLdLambda = w(j,4) + w(j,1)*x(1)*u(1) + w(j,2)*x(2)*u(2) + w(j,3)*x(3)*u(3) = 0

% Gets Lambda_bar
M1 = zeros(m+1);
M3 = zeros(m+1,1);

for i=1:m
    % each row (except end)
    M1(i,i) = 2*a(i)*u(i);
    M1(i, end) = w(timestep,i)*u(i);
    % each column (except end)
    M1(end,i) = w(timestep,i);
    M3(i,1) = -b(i)*u(i);
end
clear i
M3(end,1) = -w(timestep,4);

for i=m:-1:1
    if u(i)==0
        M1(i,:) = [];
        M1(:,i) = [];
        M3(i) = [];
    end
end

M2 = inv(M1)*M3;

step = 1;
for i=1:m
    if u(i)==0
        x_bar(i) = 0;
    else
        x_bar(i) = M2(step);
        step = step + 1;
    end
end
Lambda_bar = M2(end);

% Gets x_bar
% For largest dP add Pdiff (which is negative)
for i=1:m
    dP(1,i) = 2*a(i)*P(1,i) + b(i);
end
[~, index_max] = max(dP);
x_bar = P;
x_bar(index_max) = P(index_max) + Pdiff;



L_Js = (a(1)*x_bar(1)*x_bar(1) + b(1)*x_bar(1) + c(1))*u(1) ...
+ (a(2)*x_bar(2)*x_bar(2) + b(2)*x_bar(2) + c(2))*u(2) ...
+ (a(3)*x_bar(3)*x_bar(3) + b(3)*x_bar(3) + c(3))*u(3) ...
+ Lambda_bar*(w(timestep,4) + w(timestep,1)*x_bar(1)*u(1) + w(timestep,2)*x_bar(2)*u(2) + w(timestep,1)*x_bar(3)*u(3));
