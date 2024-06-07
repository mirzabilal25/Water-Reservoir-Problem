
clc
clear;
close all;

%% Problem Definition
  MaxIter = 200;
 
CostFunction = @(Q,Demand) waterReleaseOpFunc(Q,Demand);
 
[DemandMax,DemandMin,Evdp,Inflow] = initProblemVariabls(); 
 nVar = 120;
 VarSize=[12 10];
 VarMin = DemandMin;
 VarMax = DemandMax; 
 Storage = zeros(12,10);

 
%% Parameters of PSO
% Parameters of PSO
MaxIt = 500;
nPop =  50;
w=0.729;
%wdamp = 0.99;
c1 = 1.494; %2.05;%1.494;
c2 = 1.494; %2.05;%1.494;

%% Initialization

emptyParticle.Position = zeros(12,10);
emptyParticle.Velocity = zeros(12,10);
emptyParticle.Cost = [];
emptyParticle.Storage = zeros(12,10);
emptyParticle.Spill = zeros(12,10);

emptyParticle.Best.Position = zeros(12,10);
emptyParticle.Best.Cost = [];
emptyParticle.Best.Storage = zeros(12,10);
emptyParticle.Best.Spill = zeros(12,10);

particle = repmat(emptyParticle,nPop,1);

% Initialize GlobalBest
    GlobalBest.Cost = inf;

for i = 1:nPop
   % Generate Random Soln
   % particle(i).Position = unifrnd(VarMin, VarMax,VarSize);
   for t1 = 1: 12
       for t2 = 1:10
       %particle(i).Position(t) = unifrnd(VarMin(t), VarMax(t));
       particle(i).Position(t1,t2)=VarMin(t1,t2)+(VarMax(t1,t2)-VarMin(t1,t2)).*rand(size(VarMin(t1,t2)));
       end
   end
    
    % Intitlaize vel
    particle(i).Velocity = zeros(VarSize);
    
    %Evaluation
    [particle(i).Position,particle(i).Storage,particle(i).Spill] = CheckStorage(particle(i).Position,Evdp,Inflow,VarMax,VarMin);
    [particle(i).Cost] = CostFunction(particle(i).Position,DemandMax);
    
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

for it = 1:MaxIt
    
    for i = 1:nPop
    
        particle(i).Velocity =  w*particle(i).Velocity...
            + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position)...
            + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
        
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        %apply upper bound and lower bound limits
         %particle(i).Position = max( particle(i).Position, VarMin);
         for t1 = 1: 12
             for t2 = 1:10
                   particle(i).Position(t1,t2) = max( particle(i).Position(t1,t2), VarMin(t1,t2));
             end
         end
         for t1 = 1: 12
             for t2 = 1:10
                   particle(i).Position(t1,t2) = min( particle(i).Position(t1,t2), VarMax(t1,t2));
             end
         end
            %{
            if particle(i).Position(t) == VarMin(t)
                particle(i).Position(t)=VarMin(t)+(VarMax(t)-VarMin(t)).*rand(size(VarMin(t)));  
            end
            %}
         
         %particle(i).Position = min(particle(i).Position, VarMax);
         %for t = 1: nVar
         %   particle(i).Position(t) = min( particle(i).Position(t), VarMax(t));
         %end
        %{
         if flag == true
             for t = 1: nVar
                    particle(i).Position(t) = unifrnd(VarMin(t), VarMax(t));
             end
         end
         %}
        [particle(i).Position,particle(i).Storage,particle(i).Spill] = CheckStorage(particle(i).Position,Evdp,Inflow,VarMax,VarMin); 
        particle(i).Cost = CostFunction(particle(i).Position,DemandMax);
        
        if particle(i).Cost < particle(i).Best.Cost
           
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            particle(i).Best.Storage = particle(i).Storage;
            particle(i).Best.Spill = particle(i).Spill;
            
        end
         
    end
    %Update Global Best
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best; 
            
            end
    BestCosts(it) = GlobalBest.Cost;
    BestPositions(it) = GlobalBest;
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]); %'Storage : ' num2str(GlobalBest.Storage)
    
    %Damping Inertia coefficient
        %w = w * wdamp;
    
end
 
%% Results
figure;
semilogy(BestCosts,'LineWidth',2);
% plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Cost');
grid on;

