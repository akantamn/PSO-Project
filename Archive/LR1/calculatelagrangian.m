function [L_Js, x_bar, Lambda_bar] = calculatelagrangian(a, b, c, w, u, T, iter)
% INPUTS: u(size m), a(m), b(m), c(m), w(T,m+1), T
% OUTPUTS: L_Js % L_Js(j,:) = calculatelagrangian(a, b, c, w, u, T);

m = size(a,2); % number of variables
n = size(w,1); % number of constraints
msize = sum(u) + 1;
M1 = zeros(msize);
M3 = zeros(msize,1);


element = 1;
for i=1:m
    if u(i)==1
        M1(element,element) = 2*a(element);
        M1(element,msize) = w(element);
        M1(msize, element) = w(element);
        
        M3(element, 1) = -b(element);
        M3(msize, 1) = -w(end); % last w element
        element = element + 1;
    end
end

if det(M1)==0
    M2 = NaN;
else
    M2 = inv(M1)*M3;
end

element = 1;
for i=1:m
    if u(i)==1
        x_bar(i) = M2(element);
        element = element + 1;
    elseif u(i) ==0
        x_bar(i) = 0;
    end
end  
Lambda_bar = M2(msize);

% for i=1:m
%     if x_bar(i) < 0
%         x_bar(i) = 0; % Check xmin
%                     % Check xmax?
%     end
% end





%     L(j) = (a(1)*x(j,1)*x(j,1) + b(1)*x(j,1) + c(1))*u(j,1) ...
%     + (a(2)*x(j,2)*x(j,2) + b(2)*x(j,2) + c(2))*u(j,2) ...
%     + (a(3)*x(j,3)*x(j,3) + b(3)*x(j,3) + c(3))*u(j,3) ...
%     + Lambda(j)*(w(j,4) + w(j,1)*x(j,1)*u(j,1) + w(j,2)*x(j,2)*u(j,2)+ w(j,3)*x(j,3)*u(j,3));

%     L = (2*a(1)*x(1) + b(1) + Lambda*w(1)) * u(1) = 0
%     L = (2*a(2)*x(2) + b(2) + Lambda*w(2)) * u(2) = 0
%     L = w(3) + w(1)*x(1)*u(1) + w(2)*x(2)*u(2) = 0
    % x1, x2, Lambda

%     M1 = [2*a(1)*u(1), 0, w(1)*u(1);
%         0, 2*a(2)*u(2), w(2)*u(2);
%         w(1)*u(1), w(2)*u(2), 0];
%     M3 = [-b(1); -b(2); -w(3)];


L_Js = (a(1)*x_bar(1)*x_bar(1) + b(1)*x_bar(1) + c(1))*u(1) + (a(2)*x_bar(2)*x_bar(2) + b(2)*x_bar(2) + c(2))*u(2) + Lambda_bar*(w(3) + w(1)*x_bar(1)*u(1) + w(2)*x_bar(2)*u(2));
