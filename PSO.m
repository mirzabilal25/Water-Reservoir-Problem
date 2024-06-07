
function out = PSO(problem, params)

%% Problem Definition
 CostFunction = problem.CostFunction;
 
 nVar = problem.nVar;
 VarSize = problem.VarSize;
 VarMin = problem.VarMin;
 VarMax = problem.VarMax;
 
Demand = problem.Demand;
Evdp = problem.Evdp;
Inflow = problem.Inflow;

Lb = repmat(min(problem.Demand),1,nVar);
Ub = repmat(max(problem.Demand),1,nVar);
FEMax = problem.FEMax;
%Storage = problem.Storage;
 
%% Parameters of PSO
MaxIt = params.MaxIt;

nPop = params.nPop;
w=params.w;
wdamp = params.wdamp;
c1 = params.c1;
c2 = params.c2;
FE = 0;
FECost = [];
BestCostVal = inf;
 
%% Initialization

emptyParticle.Position = [];
emptyParticle.Velocity = [];
emptyParticle.Storage = [];
emptyParticle.Spill = [];
emptyParticle.Cost = [];
emptyParticle.Best.Position = [];
emptyParticle.Best.Cost = [];

particle = repmat(emptyParticle,nPop,1);

% Initialize GlobalBest
    GlobalBest.Cost = inf;

for i = 1:nPop
   % Generate Random Soln
   % particle(i).Position = unifrnd(VarMin, VarMax,VarSize);
   for t = 1: nVar
       %particle(i).Position(t) = unifrnd(VarMin(t), VarMax(t));
       particle(i).Position(t)=VarMin(t)+(VarMax(t)-VarMin(t)).*rand(size(VarMin(t)));
   end
    
    % Intitlaize vel
    particle(i).Velocity = zeros(VarSize);
    %particle(i).Velocity = 0.1 .* particle(i).Position;
    
    %Evaluation
    [particle(i).Position,particle(i).Storage,particle(i).Spill] = CheckStorage(particle(i).Position,Evdp,Inflow,VarMax,VarMin);
    [particle(i).Cost,FE] = CostFunction(particle(i).Position,Demand,FE);
    if (FE == 1)
        FECost(FE) = (particle(i).Cost);
    else
        FECost(FE) = (particle(i).Cost < BestCostVal)*particle(i).Cost + (particle(i).Cost >= BestCostVal)*BestCostVal;
    end
    BestCostVal = FECost(FE);
    
    %Update Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    particle(i).Best.Storage = particle(i).Storage;
    particle(i).Best.Spill = particle(i).Spill;
    
    %Update Global Best
    if particle(i).Best.Cost < GlobalBest.Cost
       GlobalBest = particle(i).Best; 
    end
    
end
% Array to hold BestCost value in each iteration
BestCosts = zeros(MaxIt,1);
%% Main Loop of PSO
it = 1;
BCost = inf;
%for it = 1:MaxIt
%while it<=MaxIt && BCost > 0
while FE <= FEMax
    
    for i = 1:nPop
    
        particle(i).Velocity =  w*particle(i).Velocity...
            + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position)...
            + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
        
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        %{
        %apply upper bound and lower bound limits
         %particle(i).Position = max( particle(i).Position, VarMin);
         for t = 1: nVar
            particle(i).Position(t) = max( particle(i).Position(t), VarMin(t));
            %{
            if particle(i).Position(t) == VarMin(t)
                particle(i).Position(t)=VarMin(t)+(VarMax(t)-VarMin(t)).*rand(size(VarMin(t)));  
            end
            %}
         end
         %particle(i).Position = min(particle(i).Position, VarMax);
         for t = 1: nVar
            particle(i).Position(t) = min( particle(i).Position(t), VarMax(t));
         end
        %{
         if flag == true
             for t = 1: nVar
                    particle(i).Position(t) = unifrnd(VarMin(t), VarMax(t));
             end
         end
         %}
        %}
        
         particle(i).Position=( particle(i).Position>Ub).*Ub+( particle(i).Position<=Ub).* particle(i).Position; 
         particle(i).Position=( particle(i).Position<Lb).*Lb+( particle(i).Position>=Lb).* particle(i).Position;
        
        
        [particle(i).Position,particle(i).Storage,particle(i).Spill] = CheckStorage(particle(i).Position,Evdp,Inflow,VarMax,VarMin); 
        [particle(i).Cost,FE] = CostFunction(particle(i).Position,Demand,FE);
        FECost(FE) = (particle(i).Cost < BestCostVal)*particle(i).Cost + (particle(i).Cost >= BestCostVal)*BestCostVal;
        BestCostVal = FECost(FE);
        
        if particle(i).Cost < particle(i).Best.Cost
           
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            particle(i).Best.Storage = particle(i).Storage;
            particle(i).Best.Spill = particle(i).Spill;
            
            %Update Global Best
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best; 
            end
            
        end
    end
    
    BestCosts(it) = GlobalBest.Cost;
    BCost = GlobalBest.Cost;
    BestPositions(it) = GlobalBest;
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]); %'Storage : ' num2str(GlobalBest.Storage)
    
    %Damping Inertia coefficient
        w = w * wdamp;
        it = it+1;
    
end
    disp(['Total Iteration: ' num2str(it-1)]);
    out.BestSolution = GlobalBest;
    out.BestCostsPSO = BestCosts;
    out.FE = FE;
    out.FECost = FECost;
    out.iter = it-1;
    %out.BestPositionsPSO = BestPositions;

end
