
function out = SA(problem, SAparams)


%% Problem Definition
CostFunction = problem.CostFunction;
 
 nVar = problem.nVar;
 VarSize = problem.VarSize;
 VarMin = problem.VarMin;
 VarMax = problem.VarMax;
 
Demand = problem.Demand;
Evdp = problem.Evdp;
Inflow = problem.Inflow;
FEMax = problem.FEMax;

%% Parameters of SA

MaxIt=SAparams.MaxIt;               % Maximum Number of Iterations
MaxSubIt=SAparams.MaxSubIt;         % Maximum Number of Sub-iterations
T0=SAparams.T0;                     % Initial Temp.
alpha=SAparams.alpha; 
Tmin =SAparams.Tmin;
FE = 0;
FECost = [];
BestCostVal = inf;
%% Initialization

% Create and Evaluate Initial Solution
%sol.Position= unifrnd(VarMin, VarMax,VarSize);
for t = 1: nVar
       sol.Position(t) = VarMin(t)+(VarMax(t)-VarMin(t)).*rand(size(VarMin(t)));
end

[sol.Position,sol.Storage,sol.Spill] = CheckStorage(sol.Position,Evdp,Inflow,VarMax,VarMin);
[sol.Cost,FE]=CostFunction(sol.Position,Demand,FE);
    if (FE == 1)
        FECost(FE) = (sol.Cost);
    else
        FECost(FE) = (sol.Cost < BestCostVal)*sol.Cost + (sol.Cost >= BestCostVal)*BestCostVal;
    end
    BestCostVal = FECost(FE);

% Initialize Best Solution Ever Found
BestSol=sol;

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

% Intialize Temp.
T=T0;

%% Main Loop of SA
it = 1;
BCost = inf;
%for it=1:MaxIt
%while it<=MaxIt && BCost > 0
while FE <= FEMax
    
    for subit=1:MaxSubIt
        
        %newsol.Position=unifrnd(VarMin, VarMax,VarSize);
        %newsol.Position=unifrnd(VarMin, VarMax,VarSize);
        for t = 1: nVar
            newsol.Position(t) = VarMin(t)+(VarMax(t)-VarMin(t)).*rand(size(VarMin(t)));
        end
        
        [newsol.Position,newsol.Storage,newsol.Spill] = CheckStorage(newsol.Position,Evdp,Inflow,VarMax,VarMin);
        [newsol.Cost,FE]=CostFunction(newsol.Position,Demand,FE);
        FECost(FE) =  (newsol.Cost < BestCostVal)*newsol.Cost + (newsol.Cost >= BestCostVal)*BestCostVal;
        BestCostVal = FECost(FE);
    
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
    BCost = BestSol.Cost;
    % Display Iteration Information
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
  
    % Update Temp.
    if T < Tmin
        break;
    else
        T=alpha*T;
    end
   it = it+1;  
end
    disp(['Total Iteration: ' num2str(it-1)]);
    out.BestSol = BestSol;
    out.BestCostsSA = BestCost;
    out.FE = FE;
    out.FECost = FECost;
    out.iter = it-1;
end
%{
%% Results
figure,
semilogy(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost SA');
grid on;

disp(['Best Sol Position SA:' num2str(BestSol.Position)]);
disp(['Solution Cost SA:' num2str(BestSol.Cost)]);
%}

