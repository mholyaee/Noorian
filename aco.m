function [BestSol,BestAnt,tau]=aco(W,TR,Ratings)
% clc;
% clear;
% close all;
gama=1.66;
%% Problem Definition

model=CreateModel(W);

CostFunction=@(tour,TR,Ratings,tau,k) TourLength(tour,model,TR,Ratings,tau,k);

nVar=model.n;


%% ACO Parameters

MaxIt=300;      % Maximum Number of Iterations

nAnt=60;        % Number of Ants (Population Size)

Q=1;

tau0=10*Q/(nVar*mean(model.D(:)));	% Initial Phromone

alpha=0.6;        % Phromone Exponential Weight
beta=0.4;         % Heuristic Exponential Weight

rho=0.02;       % Evaporation Rate


%% Initialization

eta=1./model.D;             % Heuristic Information Matrix

tau=tau0*ones(nVar,nVar);   % Phromone Matrix

BestCost=zeros(MaxIt,1);    % Array to Hold Best Cost Values

% Empty Ant
empty_ant.Tour=[];
empty_ant.Cost=[];

% Ant Colony Matrix
ant=repmat(empty_ant,nAnt,1);

% Best Ant
BestSol.Cost=inf;


%% ACO Main Loop

for it=1:MaxIt
    
    % Move Ants
    for k=1:nAnt
        
        ant(k).Tour=randi([1 nVar]);
        
        for l=2:nVar
            
            i=ant(k).Tour(end);
            %temp1=tau(i,:)
            %temp2=eta(i,:)
            myeta=eta(i,:);
            myeta(myeta>gama)=0;
            %P=tau(i,:).^alpha.*myeta.^beta;
            P=tau(i,:).^alpha.*eta(i,:).^beta;
            P(ant(k).Tour)=0;
            
            P=P/(sum(P));
            
            j=RouletteWheelSelection(P);
            
            ant(k).Tour=[ant(k).Tour j];
            
        end
        
        %ant(k).Cost=CostFunction(ant(k).Tour,model,TR,Ratings);% Why error?
        ant(k).Cost=CostFunction(ant(k).Tour,TR,Ratings,tau,k);
        
        if ant(k).Cost<BestSol.Cost
            BestSol=ant(k);
            BestAnt=k;
        end
        
    end
    
    % Update Phromones
    for k=1:nAnt
        
        tour=ant(k).Tour;
        
        %tour=[tour tour(1)]; %#ok
        
        for l=1:nVar-1
            
            i=tour(l);
            j=tour(l+1);
            
            tau(i,j)=tau(i,j)+Q/ant(k).Cost;
            
        end
        
    end
    
    % Evaporation
    tau=(1-rho)*tau;
    
    % Store Best Cost
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    % Plot Solution
%     figure(1);
%     PlotSolution(BestSol.Tour,model);
%     pause(0.01);
end

%% Results

figure;
plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

