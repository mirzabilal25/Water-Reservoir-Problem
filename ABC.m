
function out = ABC(problem, params)


%% Problem Definition

CostFunction = problem.CostFunction;
 
 nVar = problem.nVar;
 VarSize = problem.VarSize;
 VarMin = problem.VarMin;
 VarMax = problem.VarMax;
 
Demand = problem.Demand;
Evdp = problem.Evdp;
Inflow = problem.Inflow;                                      %10;               % Decision Variables Upper Bound
FEMax = problem.FEMax;

%% ABC Settings

MaxIt=params.MaxIt;                     %100;              % Maximum Number of Iterations

nPop=params.nPop;                       %100;               % Population Size (Colony Size)

nOnlooker=params.nOnlooker;             %nPop;         % Number of Onlooker Bees

L=params.L;                             %round(0.6*nVar*nPop); % Abandonment Limit Parameter (Trial Limit)

a=params.a;                             %1;                    % Acceleration Coefficient Upper Bound
FE = 0;
FECost = [];
BestCostVal = inf;
%% Initialization

% Empty Bee Structure
empty_bee.Position=[];
empty_bee.Cost=[];
empty_bee.Storage = [];
empty_bee.Spill = [];


% Initialize Population Array
pop=repmat(empty_bee,nPop,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;

% Create Initial Population
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
    
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
end

% Abandonment Counter
C=zeros(nPop,1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% ABC Main Loop
it = 1;
BCost = inf;
%for it=1:MaxIt
%while it<=MaxIt && BCost > 0
while FE <= FEMax
    % Recruited Bees
    for i=1:nPop
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee
        newbee = empty_bee;
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        
        %apply upper bound and lower bound limits
         %particle(i).Position = max( particle(i).Position, VarMin);
         for k = 1: nVar
            newbee.Position(k) = max( newbee.Position(k), VarMin(k));
         end
         %particle(i).Position = min(particle(i).Position, VarMax);
         for k = 1: nVar
            newbee.Position(k) = min( newbee.Position(k), VarMax(k));
         end
        
        % Evaluation
        [pop(i).Position,pop(i).Storage,pop(i).Spill] = CheckStorage(pop(i).Position,Evdp,Inflow,VarMax,VarMin);
        [newbee.Cost,FE]=CostFunction(newbee.Position,Demand,FE);
        FECost(FE) =  (newbee.Cost < BestCostVal)*newbee.Cost + (newbee.Cost >= BestCostVal)*BestCostVal;
        BestCostVal = FECost(FE);
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F=zeros(nPop,1);
    MeanCost = mean([pop.Cost]);
    for i=1:nPop
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
    P=F/sum(F);
    
    % Onlooker Bees
    for m=1:nOnlooker
        
        % Select Source Site
        %i=RouletteWheelSelection(P);
        r=rand;
        C=cumsum(P);
        i=find(r<=C,1,'first');

        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee
        newbee = empty_bee;
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        %apply upper bound and lower bound limits
         %particle(i).Position = max( particle(i).Position, VarMin);
         for k = 1: nVar
            newbee.Position(k) = max( newbee.Position(k), VarMin(k));
         end
         %particle(i).Position = min(particle(i).Position, VarMax);
         for k = 1: nVar
            newbee.Position(k) = min( newbee.Position(k), VarMax(k));
         end

         
        % Evaluation
        [pop(i).Position,pop(i).Storage,pop(i).Spill] = CheckStorage(pop(i).Position,Evdp,Inflow,VarMax,VarMin);
        [newbee.Cost,FE]=CostFunction(newbee.Position,Demand,FE);
        FECost(FE) =  (newbee.Cost < BestCostVal)*newbee.Cost + (newbee.Cost >= BestCostVal)*BestCostVal;
        BestCostVal = FECost(FE);
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Scout Bees
    for i=1:nPop
        if C(i)>=L
            %pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
            for k = 1: nVar
                pop(i).Position(k) = VarMin(k)+(VarMax(k)-VarMin(k)).*rand(size(VarMin(k)));
            end
            pop(i).Cost=CostFunction(pop(i).Position);
            FECost(FE) = (pop(i).Cost < BestCostVal)*pop(i).Cost + (pop(i).Cost >= BestCostVal)*BestCostVal;
            BestCostVal = FECost(FE);
            C(i)=0;
        end
    end
    
    % Update Best Solution Ever Found
    for i=1:nPop
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    BCost = BestSol.Cost;
    % Display Iteration Information
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    it = it+1;
end
    disp(['Total Iteration: ' num2str(it-1)]);
    out.BestSol = BestSol;
    out.BestCostsABC = BestCost;
    out.FE = FE;
    out.FECost = FECost;
    out.iter = it-1;
end

%{
%% Results
figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
disp(['Best Sol Position:' num2str(BestSol.Position)]);
disp(['Solution Cost:' num2str(BestSol.Cost)]);
%}