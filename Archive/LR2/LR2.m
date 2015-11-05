clear
clc

a = [0.002, 0.0025, 0.005]; % Pi^2 % P(time, unit)
b = [10, 8, 6]; % Pi
c = [500, 300, 100]; % 1
% w = [-1, -1, 5]; % P1, P2, 1
Pd = [170, 520, 1100, 330];
w(1,:) = [-1, -1, -1, Pd(1)]; % P(1), P(2), P(3), 1; Hour 1
w(2,:) = [-1, -1, -1, Pd(2)]; % Lambda*(Pload - P1 - P2 -P3); Hour 2
w(3,:) = [-1, -1, -1, Pd(3)]; % Hour 3
w(4,:) = [-1, -1, -1, Pd(4)]; % Hour 4

Pmin = [50, 100, 100];
Pmax = [600, 400, 200];
m = size(a,2); % number of variables
n = size(w,1); % number of constraints
T = size(Pd,2); % number of time periods (hours)
alpha(1) = 0.01; % if dqdLambda>0
alpha(2) = 0.002; % if dqdLambda<0

for i=1:T
    Lambda(i) = 0;
end
clear i

Lambda_next = Lambda;


for iter=1:6
    
    %% Step 1 (Primal Solution)
    % FiP Lambda, minimize for P(i), u(i)

    Lambda = Lambda_next;
    % Minimize L = (a1P1^2+c1)u1 + (a2P2^2+c2)u2 + Lambda(5-P1u1-P2u2)
        % Why is it Lambdak?
        % Minimize L = (a1P1^2+c1 - P1 Lambda)u1 + (a2P2^2+c2 - P2 Lambda)u2 + Lambda*5
        % Ignore constant term
        % Minimize L = (a1P1^2+c1 - P1 Lambda)u1 + (a2P2^2+c2 - P2 Lambda)u2
    % Minimize each term
        % Minimum is either 0 or negative (u=0 or u=1)
        % For each term
    for j=1:T
        for i=1:m % for each term (unit)
            % Get minimum at d/dPi = (2*a(i)*P(i) + b(i) + 0 - Lambda) = 0
            P(j,i) = (Lambda(j) - b(i))/(2*a(i));
            % Constrain to Pmin, Pmax
            if P(j,i)<Pmin(i)
                P(j,i) = Pmin(i);
            elseif P(j,i)>Pmax(i)
                P(j,i) = Pmax(i);
            else
                P(j,i) = P(j,i);
            end

            % Minimize cost: If cost term is negative keep it, or else use zero (u=1, u=0)
            if (a(i)*P(j,i)*P(j,i) + b(i)*P(j,i) + c(i)) + (w(j,i)*P(j,i)*Lambda(j)) > 0
                u(j,i) = 0;
                P(j,i) = 0;
            else
                u(j,i) = 1;
            end
        end
        clear i
    end
    clear j
    
    % GET L (DO WE USE IT?)
    for j=1:T
        L(j) = (a(1)*P(j,1)*P(j,1) + b(1)*P(j,1) + c(1))*u(j,1) ...
        + (a(2)*P(j,2)*P(j,2) + b(2)*P(j,2) + c(2))*u(j,2) ...
        + (a(3)*P(j,3)*P(j,3) + b(3)*P(j,3) + c(3))*u(j,3) ...
        + Lambda(j)*(w(j,4) + w(j,1)*P(j,1)*u(j,1) + w(j,2)*P(j,2)*u(j,2)+ w(j,3)*P(j,3)*u(j,3));
    end
    clear j
    
    
    % Calculate WF (Constraint) at current iteration
    Js = L;
%     for j=1:T
%         if and(and(u(j,1)==0, u(j,2)==0), u(j,3)==0)
%             Js(j) = 10000;
%         end
%     end
%     Js = sum(Js);
    for j=1:T
        WF(j) = (w(j,4) + w(j,1)*P(j,1)*u(j,1) + w(j,2)*P(j,2)*u(j,2)+ w(j,3)*P(j,3)*u(j,3));
    end
    clear j
    
    
    
    
    
    %% Step 2 (Dual Solution)
    % FiP P(i), u(i); maximize for Lambda
        % Cannot solve for maximum of q because q is unbounded wrt Lambda (?)
        % Instead, form gradient of q wrt Lambda, and adjust Lambda to move in
        % direction of increasing q
    % Get (next) Lambda
    for j=1:T
        dqdLambda(j) = (w(j,4) + w(j,1)*P(j,1)*u(j,1) + w(j,2)*P(j,2)*u(j,2)+ w(j,3)*P(j,3)*u(j,3));
        if dqdLambda(j) >= 0
            Lambda_next(j) = Lambda(j) + dqdLambda(j)*alpha(1);
        else
            Lambda_next(j) = Lambda(j) + dqdLambda(j)*alpha(2);
        end
    end
    clear j
    
    % Calculate qs at current iteration
    for j=1:T
        L(j) = (a(1)*P(j,1)*P(j,1) + b(1)*P(j,1) + c(1))*u(j,1) ...
        + (a(2)*P(j,2)*P(j,2) + b(2)*P(j,2) + c(2))*u(j,2) ...
        + (a(3)*P(j,3)*P(j,3) + b(3)*P(j,3) + c(3))*u(j,3) ...
        + Lambda(j)*(w(j,4) + w(j,1)*P(j,1)*u(j,1) + w(j,2)*P(j,2)*u(j,2)+ w(j,3)*P(j,3)*u(j,3));
    end
    clear j
    qs = L;
    qs = sum(L);
    
    %% Should be solution to Lambda, Pi from L(above)
    % Get Js (L_Js)
    for j=1:T
        
        for i=1:m
            Matrix(j,i) = 2*a(i)*P(j,i) + b(i);
        end
        
        
        
        
        PdMinusSumPU(j,:) = Pd(j) - sum(P(j,:));
        if PdMinusSumPU(j,:) < 0 
            [L_Js(j,:), P_bar(j,:), Lambda_bar(j,:)] = calculatelagrangian(a,b,c,w,u(j,:),P(j,:),PdMinusSumPU(j),j,iter); 
        elseif PdMinusSumPU(j,:) >= 0
            P_bar(j,[1:3]) = 0;
            L_Js(j,:) = 10000;
        end
    end
    clear j
    Js = L_Js;
    Js = sum(Js);
    
    
    
    %% RDG
    % qs 0->inc, Js large->dec, but qs & Js will never become equal with
    % u(i) in problem
    
    RDG = (Js-qs)/qs;
    
    %% Table
    
%     Table(iter,:) = [iter, Lambda, u, P, qs, WF, Lambda_bar, P_bar, Js, RDG];
    for j=1:T
        %Table1(j,:,iter) = [j, Lambda(j), u(j,:), P(j,:), 0, P_bar(j,:)];
        Table1(j,:,iter) = [j, Lambda(j), u(j,:), P(j,:), PdMinusSumPU(j,:), P_bar(j,:)];
    end
    clear j
    
    % Table2() = [q, Js, RDG];
    Table2(iter,:) = [qs, Js, RDG];

clear PdMinusSumPU
end


% Need to fiP Js, RDG

disp([' iter, Lambda, u(1), u(2), u(3), P(1), P(2), P(3), WF, P_bar(1), P_bar(2), P_bar(3)']);
Table1
disp(['Table2 q (Primal Solution), Js (Dual Solution), RDG (Relative Duality Gap)']);
Table2
