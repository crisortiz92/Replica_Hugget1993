%                       UNIVERSIDAD SAN FRANCISCO DE QUITO
%               
%                                     PART II
%
%                     REPLICATION OF FINDINGS BY HUGGET (1993)
%                                
%           WE REPLICATE BOTH TABLES AND GRAPHS OF THE PAPER : THE
%           RISK-FREE RATE IN HETEROGENEOUS-AGENT INCOMPLETE-INSURANCE
%           ECONOMIES. WE USE ENDOGENOUS GRID METHOD FOR OPTIMIZATION
%                           
%           OPTIMIZATION: USING MONTECARLO SIMULATION
%
%                                CRISTHIAN ORTIZ
%                                PAUL PONCE
%                                KEVIN ROJAS
%                                SANTIAGO SANDOVAL
%       
%

% tidying up
clear all
close all
clc
tic
disp('Huggett model; credit market clearing');


%  ECONOMIC PARAMETERS

sigma  = 1.5;                                           % risk aversion              
nop = 6;                                                % number of model periods in a year,1 means annual, 4 means quarterly model
betaAnnual   = 0.96;                                    % subjective discount factor annual
beta = betaAnnual^(1/nop);                              % subjective discount factor in model
pee = 0.925;                                            % probability of going from Employed to Employed
puu = 0.5;                                              % probability of goinge from Unemployed to Unemployed
prob   = [pee (1-pee);(1-puu) puu];                     % prob(i,j) = probability (s(t+1)=sj | s(t) = si)
ns = length(prob);                                      % number of employment / productivity states
theta  = 0.1;                                           % non-interest income if unemployed
wage   = 1.00;                                          % non-interest income if employed
income = [wage, theta*wage];                            % stack it into vector
Rh = 1/beta-0.001;                                      % initial high gross interest rate 
phiadhoc = -2;                                          % ad hoc borrowing constraint parameter 
flagBC=1;                                               % =0 natural; =1 ad hoc
%
%  TECHNICAL PARAMETERS
  
maxa = 35;                                              % maximum value of asset grid   high value necessary
na = 100;                                               % # asset nodes
flaggrida = 2;                                          % 1=linear grid, 2=exponential grid, 3= double exponential
maxiterR   = 10;                                        % maximum number of interest rate iterations
tolEE   = 1.e-5;                                        % tolerance level that measures convergence of EE
maxiterEE = 5000;                                       % maximum number of EE iterations 
naF=500;                                                % # assets for optimization
flaggridaF=2;                                           % Type of grid for optimization
tolpop = 1.e-10;                                        % tolerance level pop density
maxiterpop = 10000;                                     % max iteration for matrix transition
pop=ones(1,naF*ns)/(naF*ns);                            % initial population, uniform, anything with mass 1 will do.
tolED = 1.e-3;                                          % tolearance level for excess demand in credit market

maxiterED=30;                                           % max iteration for credit market equilibrium
numberHome=1000;
numberTime=500;
MaxTime = 10000;
tolMont=0.001;


R = Rh;     % Start with ih
critED=1;
iterED=1;
    
    phi=phiadhoc;
  
   
    % The grid for savings that are chosen
    mina=0.001+phi;         % minimum value of asset grid , depends on phi
    grida=makegrid(flaggrida,na,mina,maxa)';    % asset (tommorow) grid
    gridaM=grida*ones(1,ns);        % in matrix form (needed later)

    coh_grid=R*grida*ones(1,ns);    % everyone gets interest income
    coh_grid=coh_grid+ones(na,1)*income;    % labor income according to state

    cpr=coh_grid-gridaM;    % guess for tomorrow's consumption
    cpr(cpr<.01*min(income))=.01*min(income);  % avoid negative c

    
    %%%%%        SOLVING VALUE & POLICY FUNCTIONS
    
    critEE=1; % initialize crit
    iterEE=0;   % start loop counter
    while ((critEE>tolEE) && (iterEE<maxiterEE))
        cprO=cpr;
        iterEE=iterEE+1;
        EMUpr=cpr.^(-sigma)*prob'; % 1. step: Compute expected marginal utility prime (tomorrow)
        EMUpr=beta*R*EMUpr;     % 2. discount it
        c=EMUpr.^(-1/sigma);    % 3. invert Euler eqn
        coh=c+gridaM;           % 4. Cash on hand implied
        for j=1:ns  % interpolate consumption back onto our fixed coh_grid
            cpr(:,j)=interp1(coh(:,j),c(:,j),coh_grid(:,j),'pchip','extrap');
        end
        if flagBC==1    %check BC
            cpr=min(cpr,coh_grid-phi);
        end
        critEE=max(max(abs(cprO-cpr)./(1+cpr)));
    end
    
    
   %%%%%%%   COMPUTING MONTECARLO
    
   sav=coh_grid-cpr;                                                %Policy Function
   rng(111);                                                        %Fixing seed
   assets=phiadhoc + (maxa - phiadhoc) * rand(numberHome,1);        % First vector of assets
   endowments= theta + (wage - theta).* (randi(2,numberHome,1)-1);  % First vector of endowments
   MatEndowments= zeros(numberTime,numberHome);                     % Matrix for endowments
   MC = dtmc(prob);                                                 % Transition as an object
   Matassets = zeros(MaxTime,numberHome);                           % Create array for assets
   Matassets(1,:) = assets;                                         % Assigning first element
  
   statenames = [0.1  1];
   for i=1:numberHome                                               %Creating Markov Chains for endowments
        X = simulate(MC,numberTime-1);
        MatEndowments(:,i) = X';   
   end
   
   itMont = 1;                                                      %Iteration of convergence
   convMont = 1;                                                    %initial for convergence
   
   while(convMont > tolMont)
   
      for j=1:numberHome
           Matassets(itMont+1,j) =  interp1(grida,sav(:,MatEndowments(itMont,j)),Matassets(itMont,j),'pchip','extrap');
      end
  
   % Measures
   
   Konemean = mean(Matassets(itMont,:));                            % mean for a present
   Konesd = std(Matassets(itMont,:));                               % sd for a present
   Eonemean = std(MatEndowments(itMont,:));                         % mean for endowments present
   Eonestd = std(MatEndowments(itMont,:));                          % sd for endowments present
   Ktwomean = mean(Matassets(itMont+1,:));                          % mean for apri future
   Ktwosd = std(Matassets(itMont+1,:));                              % sd for endowments 
   Etwomean = mean(MatEndowments(itMont+1,:));                      % mean for endowments future
   Etwosd = std(MatEndowments(itMont+1,:));                         % sd for endowments future
   
   itMont = itMont + 1;                                             % count of iterations
   
   convMont = abs(Ktwomean-Konemean);                               % criteria for convergence of mean a
   convEndow = abs(Etwomean-Eonemean);                              % criteria for convergence of mean endow
   convsdMont = abs(Ktwosd-Konesd);                                 % criteria for convergence of sd a
   convsdEndow = abs(Etwosd-Eonestd);                                % criteria for convergence of sd endow
   
   end
  
   
   % Computing consume
   
   Cgen = zeros(1, numberHome);                                     % consume for last period using interpol
   
   for i=1:numberHome
      Cgen(1,i) = interp1(grida,c(:,MatEndowments(itMont-1,i)),Matassets(itMont-1,i),'pchip','extrap'); 
   end
   
   % Computing q
   
   qVect = zeros(1, numberHome);                                    % calculating bondprices
   
   for i=1:numberHome
      qVect(1,i) = (Matassets(itMont-1, i) + MatEndowments(itMont-1, i) - Cgen(1,i)) / Matassets(itMont, i); 
   end
   
   % Estimating via mean
   q = mean(qVect);
   disp('The bond price equilibrium for model periods is:')
   disp(q)
   
   R = 1/q;
   disp('The credit market equilibrium at annualized interest rate is:')
   disp((R^6-1)*100)
   
  toc
  
AssetsEnd = Matassets(itMont,:);
EndowmentsEnd = MatEndowments(itMont,:);
EndowmentsInit = MatEndowments(1,:);
MatPlots = [EndowmentsInit; EndowmentsEnd; AssetsEnd];
Graficos = MatPlots';
AssetsEndOrder = sort(Graficos);
figure
hold on
plot(sort(Graficos(1:498,3)), 'r:','LineWidth',1.);
plot(sort(Graficos(499:1000,3)), 'b:','LineWidth',1.);
xlim([1 10])
ylim([-2 4])
title('Fig. 1 Optimal decision rule.' )
xlabel('Credit Level = a')
legend('SAVINGS EMPLOYED','SAVINGS UNEMPLOYED' ,'Location','SouthEast')
