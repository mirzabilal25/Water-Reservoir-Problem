clc
clear;
close all;

%% Problem Definition
 MaxIter = 80;
 popSize = 15;
 CostFunction = @(Q,Demand) waterReleaseOpFunc(Q,Demand);
 
[Demand,Evdp,Inflow] = initProblemVariabls(); 
 nVar = 12;
 VarSize=[1 nVar];
 VarMin = zeros(1,nVar);
 VarMax = Demand; 
 Storage = zeros(1,nVar);

%% Parameters of SA

MaxIt=MaxIter;        % Maximum Number of Iterations
MaxSubIt= 65;       % Maximum Number of Sub-iterations
T0=0.9;          % Initial Temp.
alpha=0.8;
Tmin = 1e-70;
c1 = 0.91;

%% Initialization

% Create and Evaluate Initial Solution
%sol.Position= unifrnd(VarMin, VarMax,VarSize);
for t = 1: nVar
       sol.Position(t) = VarMin(t)+(VarMax(t)-VarMin(t)).*rand(size(VarMin(t)));
end

[sol.Position,sol.Storage,sol.Spill] = CheckStorage(sol.Position,Evdp,Inflow,VarMax,VarMin);
sol.Cost=CostFunction(sol.Position,Demand);

% Initialize Best Solution Ever Found
BestSol=sol;

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

% Intialize Temp.
T=T0;

%% Main Loop of SA

for it=1:MaxIt
    
    for subit=1:MaxSubIt
        
        %newsol.Position=unifrnd(VarMin, VarMax,VarSize);
        for t = 1: nVar
            newsol.Position(t) = VarMin(t)+(VarMax(t)-VarMin(t)).*rand(size(VarMin(t)));
        end
        
        [newsol.Position,newsol.Storage,newsol.Spill] = CheckStorage(newsol.Position,Evdp,Inflow,VarMax,VarMin);
        newsol.Cost=CostFunction(newsol.Position,Demand);
    
        if newsol.Cost<=sol.Cost % If NEWSOL is better than SOL
            sol=newsol;
            %break;
        else
                DELTA=(newsol.Cost-sol.Cost);
                if  (DELTA < 1e-6)                                    % If NEWSOL is NOT better than SOL
                    sol=newsol;
            
                else
                    P=exp(DELTA/T);
                    if rand < P
                        sol=newsol;
                    end
            
                end
        end
        
        % Update Best Solution Ever Found
        if sol.Cost <= BestSol.Cost
            BestSol = sol;
        end
        
    end     
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
  
    % Update Temp.
    if T < Tmin
        break;
    else
        T=alpha*T;
    end
     
end

%% Results
figure,
semilogy(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost SA');
grid on;

disp(['Best Sol Position SA:' num2str(BestSol.Position)]);
disp(['Solution Cost SA:' num2str(BestSol.Cost)]);

