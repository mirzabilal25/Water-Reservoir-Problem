
function out = GA(problem, GAparams)


%% Problem Definition

CostFunction = problem.CostFunction;
 
 nVar = problem.nVar;
 VarSize = problem.VarSize;
 VarMin = problem.VarMin;
 VarMax = problem.VarMax;
 
Demand = problem.Demand;
Evdp = problem.Evdp;
Inflow = problem.Inflow;                            %10;                % Upper Bound of Variables
FEMax = problem.FEMax;

%% GA Parameters

MaxIt=GAparams.MaxIt;           %300;     % Maximum Number of Iterations

nPop=GAparams.nPop;             %100;       % Population Size

pc=GAparams.pc;                 %0.7;                 % Crossover Percentage
nc=GAparams.nc;                 %2*round(pc*nPop/2);  % Number of Offsprings (also Parnets)
gamma=GAparams.gamma;           %0.4;              % Extra Range Factor for Crossover

pm=GAparams.pm;                 %0.3;                 % Mutation Percentage
nm=GAparams.nm;                 %round(pm*nPop);      % Number of Mutants
mu=GAparams.mu;                 %0.1;                 % Mutation Rate

% RouletteWheelSelection
beta=GAparams.beta;             %8;  % Selection Pressure

FE = 0;
FECost = [];
BestCostVal = inf;
%pause(0.01);

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Storage = [];
empty_individual.Spill = [];
pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
    % Initialize Position
    %pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    for k = 1: nVar
       pop(i).Position(k) = VarMin(k)+(VarMax(k)-VarMin(k)).*rand(size(VarMin(k)));
    end
    
    % Evaluation
    [pop(i).Position,pop(i).Storage,pop(i).Spill] = CheckStorage(pop(i).Position,Evdp,Inflow,VarMax,VarMin);
    [pop(i).Cost,FE]=CostFunction(pop(i).Position,Demand,FE);
    if (FE == 1)
        FECost(FE) = (pop(i).Cost);
    else
        FECost(FE) = (pop(i).Cost < BestCostVal)*pop(i).Cost + (pop(i).Cost >= BestCostVal)*BestCostVal;
    end
    BestCostVal = FECost(FE);
    
end

% Sort Population
Costs=[pop.Cost];
[Costs, SortOrder]=sort(Costs);
pop=pop(SortOrder);

% Store Best Solution
BestSol=pop(1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

% Store Cost
WorstCost=pop(end).Cost;

%% Main Loop
it = 1;
BCost = inf;
%for it=1:MaxIt
%while it<=MaxIt && BCost > 0
while FE <= FEMax
    % Calculate Selection Probabilities
     P=exp(-beta*Costs/WorstCost);
     P=P/sum(P);
     
     % Crossover
    popc=repmat(empty_individual,nc/2,2);
    for k=1:nc/2
        i1=RouletteWheelSelection(P);
        i2=RouletteWheelSelection(P);
        
        p1=pop(i1);
        p2=pop(i2);
        
        [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position,gamma,VarMin,VarMax);
        
        [popc(k,1).Position,popc(k,1).Storage,popc(k,1).Spill] = CheckStorage(popc(k,1).Position,Evdp,Inflow,VarMax,VarMin);
        [popc(k,1).Cost,FE]=CostFunction(popc(k,1).Position,Demand,FE);
        FECost(FE) = (popc(k,1).Cost < BestCostVal)*popc(k,1).Cost + (popc(k,1).Cost >= BestCostVal)*BestCostVal;
        BestCostVal = FECost(FE);
        
        
        [popc(k,2).Position,popc(k,2).Storage,popc(k,2).Spill] = CheckStorage(popc(k,2).Position,Evdp,Inflow,VarMax,VarMin);
        [popc(k,2).Cost,FE]=CostFunction(popc(k,2).Position,Demand,FE);
        FECost(FE) = (popc(k,2).Cost < BestCostVal)*popc(k,2).Cost + (popc(k,2).Cost >= BestCostVal)*BestCostVal;
        BestCostVal = FECost(FE);
        
    end
    popc=popc(:);
    
    % Mutation
    popm=repmat(empty_individual,nm,1);
    for k=1:nm
        
        % Select Parent
        i=randi([1 nPop]);
        p=pop(i);
        
        % Apply Mutation
        popm(k).Position=Mutate(p.Position,mu,VarMin,VarMax);
        
        % Evaluate Mutant
        [popm(k).Position,popm(k).Storage,popm(k).Spill] = CheckStorage(popm(k).Position,Evdp,Inflow,VarMax,VarMin);
        [popm(k).Cost,FE]=CostFunction(popm(k).Position,Demand,FE);
        FECost(FE) = (popm(k).Cost < BestCostVal)*popm(k).Cost + (popm(k).Cost >= BestCostVal)*BestCostVal;
        BestCostVal = FECost(FE);
        
    end
    
     % Create Merged Population
     pop=[pop
         popc
         popm];
    
    % Sort Population
    Costs=[pop.Cost];
    [Costs, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
    
    % Update Worst Cost
    WorstCost=max(WorstCost,pop(end).Cost);
    
    % Truncation
    pop=pop(1:nPop);
    Costs=Costs(1:nPop);
    
    % Store Best Solution Ever Found
    BestSol=pop(1);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    BCost = BestSol.Cost;
    % Show Iteration Information
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    it = it+1;
end
    disp(['Total Iteration: ' num2str(it-1)]);
    out.BestSol = BestSol;
    out.BestCostsGA = BestCost;
    out.FE = FE;
    out.FECost = FECost;
    out.iter = it-1;
end
%{
%% Results
figure;
semilogy(BestCost,'LineWidth',2);
% plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Cost');
grid on;
%}