%                       UNIVERSIDAD SAN FRANCISCO DE QUITO
%               
%                                     PART I
%
%                     REPLICATION OF FINDINGS BY HUGGET (1993)
%                                
%           WE REPLICATE BOTH TABLES AND GRAPHS OF THE PAPER : THE
%           RISK-FREE RATE IN HETEROGENEOUS-AGENT INCOMPLETE-INSURANCE
%           ECONOMIES. WE USE ENDOGENOUS GRID METHOD FOR OPTIMIZATION
%                           
%           OPTIMIZATION: USING DENSITY FUNCTION AND BISECTION METHOD
%
%                                CRISTIAN ORTIZ
%                                PAUL PONCE
%                                KEVIN ROJAS
%                                SANTIAGO SANDOVAL
%       
%


% Formal beggining of the document
clear all
close all
clc

tic 
disp('Huggett model; credit market clearing');


%%%%%%%  ECONOMIC PARAMETERS

sigma  = 1.5;                               % risk aversion              
nop = 6;                                    % number of model periods in a year,1 means annual, 4 means quarterly model
betaAnnual   = 0.96;                        % subjective discount factor annual
beta = betaAnnual^(1/nop);                  % subjective discount factor in model
pee = 0.925;                                % probability of going from Employed to Employed
puu = 0.5;                                  % probability of goinge from Unemployed to Unemployed
prob   = [pee (1-pee);(1-puu) puu];         % prob(i,j) = probability (s(t+1)=sj | s(t) = si)
ns = length(prob);                          % number of employment / productivity states
theta  = 0.1;                               % non-interest income if unemployed
wage   = 1.00;                              % non-interest income if employed
income = [wage, theta*wage];                % stack it into vector
Rh = 1/beta-0.001;                          % initial high gross interest rate 
phiadhoc = -2;                              % ad hoc borrowing constraint parameter
flagBC = 1;                                 

%%%%%%%  TECHNICAL PARAMETERS
  
maxa = 35;                                  % maximum value of asset grid   high value necessary
na = 100;                                   % # asset nodes
flaggrida = 2;                              % 1=linear grid, 2=exponential grid, 3= double exponential
maxiterR   = 10;                            % maximum number of interest rate iterations
tolEE   = 1.e-5;                            % tolerance level that measures convergence of EE
maxiterEE = 5000;                           % maximum number of EE iterations 
naF=500;                                    % # assets for optimization
flaggridaF=2;                               % Type of grid for optimization
tolpop = 1.e-10;                            % tolerance level pop density
maxiterpop = 10000;                         % max iteration for matrix transition
pop=ones(1,naF*ns)/(naF*ns);                % initial population, uniform, anything with mass 1 will do.
tolED = 1.e-3;                              % tolearance level for excess demand in credit market
flaghunt=1;                                 % flag for convergence of bisection
maxiterED=30;                               % max iteration for credit market equilibrium


R = Rh;     % Start of interest rate (from theory)
critED=1;   % Compared to tolerance
iterED=0;   % Number of iterations for update R

while ((critED>tolED) && (iterED<maxiterED))
    iterED=iterED+1;
    if flagBC==1
        phi=phiadhoc;             % Rename
    else
        phi=-min(income)/(R-1);   % natural borrowing limit
    end
    
    % The grid for savings that are chosen
    mina=0.001+phi;                             % minimum value of asset grid , depends on phi
    grida=makegrid(flaggrida,na,mina,maxa)';    % asset (tommorow) grid
    gridaM=grida*ones(1,ns);                    % in matrix form (needed later)

    coh_grid=R*grida*ones(1,ns);                % everyone gets interest income
    coh_grid=coh_grid+ones(na,1)*income;        % labor income according to state

    cpr=coh_grid-gridaM;                        % guess for tomorrow's consumption
    cpr(cpr<.01*min(income))=.01*min(income);   % avoid negative c

    
    
    %%%%%%%%        SOLVING VALUE & POLICY FUNCTIONS
    
    critEE=1; % initialize crit
    iterEE=0;   % start loop counter
    while ((critEE>tolEE) && (iterEE<maxiterEE))
        cprO=cpr;
        iterEE=iterEE+1;
        EMUpr=cpr.^(-sigma)*prob';                  % 1. step: Compute expected marginal utility prime (tomorrow)
        EMUpr=beta*R*EMUpr;                         % 2. discount it
        c=EMUpr.^(-1/sigma);                        % 3. invert Euler 
        coh=c+gridaM;                               % 4. Cash on hand implied
        
        for j=1:ns                                  % interpolate consumption back onto our fixed coh_grid
            cpr(:,j)=interp1(coh(:,j),c(:,j),coh_grid(:,j),'pchip','extrap');
        end
        if flagBC==1    %check BC
            cpr=min(cpr,coh_grid-phi);
        end 
        critEE=max(max(abs(cprO-cpr)./(1+cpr)));    % convergence criteria
    end
    
    
    %%%%%%%%   COMPUTING THE STEADY STATE POPULATION DENSITY 
    

    sav=coh_grid-cpr;                               %assigning savings
    gridaF=makegrid(flaggridaF,naF,mina,maxa)';     % fine asset grid
    savF=zeros(naF,ns);
    for j=1:ns                                      % interpolate savings function onto fine grid
        savF(:,j)=interp1(grida,sav(:,j),gridaF,'pchip','extrap');
    end
    savI=zeros(naF,ns);                             % find the index of gridaF that corresponds to this savings choice
    for j=1:ns
        for ia=1:naF
            [temp, savI(ia,j)]=min(abs(savF(ia,j)-gridaF));
        end
    end
    % Create the (sparse) transition matrix
    
    rowI=zeros(naF*ns*ns,1);                        % row index
    colI=zeros(naF*ns*ns,1);                        % col index
    valI=zeros(naF*ns*ns,1);                        % matrix value
    ind=0;
    for j=1:ns                                      % states today
        for ia=1:naF                                % asset today
            for jj=1:ns                             % state tomorrow
                ind=ind+1;
                rowI(ind)=(j-1)*naF+ia;
                colI(ind)=(jj-1)*naF+savI(ia,j);
                valI(ind)=prob(j,jj);
            end
        end
    end
    trans=sparse(rowI,colI,valI,naF*ns,naF*ns);     % need to tell matlab length of rows and cols 
    
    % 6. Solve invariant density
    
    critpop=1;                      
    iter=1;
    while ((critpop>tolpop) && (iter<maxiterpop))
        popnew=pop*trans;                           % apply transition
        critpop=max(abs(popnew-pop));
        pop=popnew;
        iter=iter+1;
    end
    
    %%%%%%                BISECTION METHOD BASED ON EQUATIONS
    meanA=pop*savF(:);
    if flaghunt==1
        if meanA>0                                      % still not bracketed true R
            Rh=Rh-0.01;
            R=Rh;
        else                                            % brackted true R for first time
            Rl=R;
            Rh=Rh+0.01;
            R=(Rh+Rl)/2;
            flaghunt=0;
        end
    else                                                % bisect
        if meanA>0
            Rh=R;
        else
            Rl=R;
        end
        R=(Rh+Rl)/2;
    end
    disp('new interest rate')
    disp(R)
    critED=abs(meanA);        
end

bondprice=1/R;
disp('The bond price equilibrium for model periods is:')              %Display bond price
disp(bondprice)
    

disp('The credit market equilibrium at annualized interest rate is:') %Display equilibrium
disp((R^6-1)*100)



% SOME ADDITIONAL RESULTS FOR GRAPHS
popSS=reshape(pop',naF,ns);                                           % density: col1 E ; col2 U
popCUM=cumsum(popSS);
indRich=find(popSS(:,1)>0,1,'last');                                  % identifies the richest agent
ub=indRich+100;

%Plots
% (GROSS) SAVINGS
figure
hold on
plot(gridaF(1:ub),savF(1:ub,1),'b-',gridaF(1:ub),gridaF(1:ub),'r:','LineWidth',1.)
plot(gridaF(1:ub),savF(1:ub,2),'b-',gridaF(1:ub),gridaF(1:ub),'r:','LineWidth',1.)
xlim([-2 4])
ylim([-2 4])
title('Fig. 1 Optimal decision rule.')
xlabel('Credit Level = a')
legend('SAVINGS EMPLOYED','45 DEGREE LINE','SAVINGS UNEMPLOYED','Location','SouthEast')

% WEALTH DISTRIBUTION

figure
plot(gridaF(1:ub),popCUM(1:ub,1),gridaF(1:ub),popCUM(1:ub,2),'r:','LineWidth',1.0)
xlim([-2 0.97])
ylim([0 1])
xlabel('Credit Level = a')
title('Fig.2 Stationary distribution')
legend('CDF EMPLOYED','CDF UNEMPLOYED','Location','NorthWest')

toc
