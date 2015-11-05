clear
clc

a = [0.25, 0.255]; % xi^2
b = [0, 0]; % xi
c = [15, 15]; % 1
w = [-1, -1, 5]; % x1, x2, 1
xmin = [0, 0];
xmax = [10, 10];
m = size(a,2); % number of variables
n = size(w,1); % number of constraints
alpha(1) = 0.2; % if dqdLambda>0
alpha(2) = 0.005; % if dqdLambda<0

% L = F + Lambda(W)

% Get Primal Solution
    % min(f)
    % q = min L
% Get Dual Solution
    % max{Lambda} min{x} (f + Sum Lambdai hi)
    % qs = max q

Lambda = 0;
Lambda_next = Lambda;

disp([' iter, Lambda, u(1), u(2), x(1), x(2), qs, WF, Lambda_bar, x_bar(1), x_bar(2), Js, RDG']);

for iter=1:8
    
    %% Step 1 (Primal Solution)
    % Fix Lambda, minimize for x(i), u(i)

    Lambda = Lambda_next;
    % Minimize L = (a1x1^2+c1)u1 + (a2x2^2+c2)u2 + Lambda(5-x1u1-x2u2)
        % Why is it Lambdak?
        % Minimize L = (a1x1^2+c1 - x1 Lambda)u1 + (a2x2^2+c2 - x2 Lambda)u2 + Lambda*5
        % Ignore constant term
        % Minimize L = (a1x1^2+c1 - x1 Lambda)u1 + (a2x2^2+c2 - x2 Lambda)u2

    % Minimize each term
        % Minimum is either 0 or negative (u=0 or u=1)
        % For each term
    for i=1:m % for each term
        % Get minimum at d/dxi = (2*a(i)*x(i) + b(i) + 0 - Lambda) = 0
        x(i) = (Lambda - b(i))/(2*a(i));
        % Constrain to xmin, xmax
        if x(i)<xmin(i)
            x(i) = xmin(i);
        elseif x(i)>xmax(i)
            x(i) = xmax(i);
        else
            x(i) = x(i);
        end

        % if cost term is negative keep it, or else use zero (u=1, u=0)
        if (a(i)*x(i)*x(i) + b(i)*x(i) + c(i) + w(i)*x(i)*Lambda) > 0
            u(i) = 0;
        else
            u(i) = 1;
        end
    end
    
    L = (a(1)*x(1)*x(1) + b(1)*x(1) + c(1))*u(1) + (a(2)*x(2)*x(2) + b(2)*x(2) + c(2))*u(2) + Lambda*(w(3) + w(1)*x(1)*u(1) + w(2)*x(2)*u(2));
    Js = L;
    if and(u(1)==0, u(2)==0)
        Js = 50;
    end
    WF = (w(3) + w(1)*x(1)*u(1) + w(2)*x(2)*u(2));

    %% Step 2 (Dual Solution)
    % Fix x(i), u(i); maximize for Lambda
        % Cannot solve for maximum of q because q is unbounded wrt Lambda (?)
        % Instead, form gradient of q wrt Lambda, and adjust Lambda to move in
        % direction of increasing q
    dqdLambda = (w(3) + w(1)*x(1)*u(1) + w(2)*x(2)*u(2));
    if dqdLambda >= 0
        Lambda_next = Lambda + dqdLambda*alpha(1);
    else
        Lambda_next = Lambda + dqdLambda*alpha(2);
    end
    
    L = (a(1)*x(1)*x(1) + b(1)*x(1) + c(1))*u(1) + (a(2)*x(2)*x(2) + b(2)*x(2) + c(2))*u(2) + Lambda*(w(3) + w(1)*x(1)*u(1) + w(2)*x(2)*u(2));
    qs = L;
    
    %% Should be solution to Lambda, xi from L(above)
        % J=L
        % Js = minL
        % Use Lambda(prev) or Lambda_next ???
        % Keep ui as solved now, solve xi, Lambda for min L
    % u(i) known
% %     for i=1:m % Solve minimum for each cost term
% %         % Use Lambda_next
% %         x_bar(i) = (Lambda_next - b(i))/(2*a(i));
% %     end
% %     for i=1:n % Solve minimum for each constraint term
% %         
% %     end

%     L = (2*a(1)*x(1) + b(1) + Lambda*w(1)) * u(1) = 0
%     L = (2*a(2)*x(2) + b(2) + Lambda*w(2)) * u(2) = 0
%     L = w(3) + w(1)*x(1)*u(1) + w(2)*x(2)*u(2) = 0
    % x1, x2, Lambda

%     M1 = [2*a(1)*u(1), 0, w(1)*u(1);
%         0, 2*a(2)*u(2), w(2)*u(2);
%         w(1)*u(1), w(2)*u(2), 0];
%     M3 = [-b(1); -b(2); -w(3)];
    
    
    
%     if and(u(1)==0, u(2)==0)
%         Lambda_bar = 0;
%         x_bar(1) = 0;
%         x_bar(2) = 0;
%     elseif u(1)==0
%         x_bar(1) = 0;
%         % x2, Lambda
%         M1 = [2*a(2)*u(2), w(2)*u(2);
%         w(2)*u(2), 0];
%         M3 = [-b(2); -w(3)];
%         M2 = inv(M1)*M3;
%         x_bar(1) = 0;
%         x_bar(2) = M2(1);
%         Lambda_bar = M2(2);
%     elseif u(2)==0
%         x_bar(2) = 0;
%         % x1, Lambda
%         M1 = [2*a(1)*u(1), w(1)*u(1);
%         w(1)*u(1), 0];
%         M3 = [-b(1); -w(3)];
%         M2 = inv(M1)*M3;
%         x_bar(1) = M2(1);
%         x_bar(2) = 0;
%         Lambda_bar = M2(2);
%     elseif and(u(1)==1, u(2)==1)
%         M1 = [2*a(1)*u(1), 0, w(1)*u(1);
%         0, 2*a(2)*u(2), w(2)*u(2);
%         w(1)*u(1), w(2)*u(2), 0];
%         M3 = [-b(1); -b(2); -w(3)];
%         M2 = inv(M1)*M3;
%         x_bar(1) = M2(1);
%         x_bar(2) = M2(2);
%         Lambda_bar = M2(3);
%     end
    
    %% Js
%     L_Js = (a(1)*x_bar(1)*x_bar(1) + b(1)*x_bar(1) + c(1))*u(1) + (a(2)*x_bar(2)*x_bar(2) + b(2)*x_bar(2) + c(2))*u(2) + Lambda_bar*(w(3) + w(1)*x_bar(1)*u(1) + w(2)*x_bar(2)*u(2));
    [L_Js, x_bar, Lambda_bar] = calculatelagrangian(a,b,c,w,u,1);
    Js = L_Js;
    if and(u(1)==0, u(2)==0)
        Js = 50;
    end
    
    %% RDG
    % qs 0->inc, Js large->dec, but qs & Js will never become equal with
    % u(i) in problem
    
    RDG = (Js-qs)/qs;
    
    %% Table
    
    Table(iter,:) = [iter, Lambda, u, x, qs, WF, Lambda_bar, x_bar, Js, RDG];
    
end

Table
