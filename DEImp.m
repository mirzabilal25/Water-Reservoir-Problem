function outDE = DEImp(problem, DEparams) 
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

%% DE Parameters

MaxIt=DEparams.MaxIt;      % Maximum Number of Iterations

nPop=DEparams.nPop;        % Population Size

beta_min=DEparams.beta_min;   % Lower Bound of Scaling Factor
beta_max=DEparams.beta_max;   % Upper Bound of Scaling Factor

pCR=DEparams.pCR;        % Crossover Probability
FE = 0;
FECost = [];
BestCostVal = inf;
%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Storage = [];
empty_individual.Spill = [];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

for i=1:nPop

    %pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    for k = 1: nVar
       pop(i).Position(k) = VarMin(k)+(VarMax(k)-VarMin(k)).*rand(size(VarMin(k)));
    end
    [pop(i).Position,pop(i).Storage,pop(i).Spill] = CheckStorage(pop(i).Position,Evdp,Inflow,VarMax,VarMin);
    [pop(i).Cost,FE]=CostFunction(pop(i).Position,Demand,FE);
    if (FE == 1)
        FECost(FE) = (pop(i).Cost);
    else
        FECost(FE) = (pop(i).Cost < BestCostVal)*pop(i).Cost + (pop(i).Cost >= BestCostVal)*BestCostVal;
    end
    BestCostVal = FECost(FE);
    
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
    
end

BestCost=zeros(MaxIt,1);

%% DE Main Loop
it = 1;
BCost = inf;
%for it=1:MaxIt
%while it<=MaxIt && BCost > 0
while FE <= FEMax
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
        %beta=unifrnd(beta_min,beta_max);
        %beta=unifrnd(beta_min,beta_max,VarSize);
        beta=beta_min+(beta_max-beta_min).*rand(size(beta_min));
        
        y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
        
        for k = 1: nVar
            y(k) = max( y(k), VarMin(k));
         end
         %particle(i).Position = min(particle(i).Position, VarMax);
         for k = 1: nVar
            y(k) = min( y(k), VarMax(k));
         end
        
        %y = max(y, VarMin);
		%y = min(y, VarMax);
		
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        NewSol = empty_individual;
        NewSol.Position=z;
        [NewSol.Position,NewSol.Storage,NewSol.Spill] = CheckStorage(NewSol.Position,Evdp,Inflow,VarMax,VarMin);
        [NewSol.Cost,FE]=CostFunction(NewSol.Position,Demand,FE);
        FECost(FE) =  (NewSol.Cost < BestCostVal)*NewSol.Cost + (NewSol.Cost >= BestCostVal)*BestCostVal;
        BestCostVal = FECost(FE);
        
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
            end
        end
        
    end
    
    % Update Best Cost
    BestCost(it)=BestSol.Cost;
    BCost = BestSol.Cost;
    % Show Iteration Information
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]); %'Storage : ' num2str(BestSol.Storage)
    it = it+1;
    
end
    disp(['Total Iteration: ' num2str(it-1)]);
    outDE.BestSol = BestSol;
    outDE.BestCostsDE = BestCost;
    outDE.FE = FE;
    outDE.FECost = FECost;
    outDE.iter = it-1;
end