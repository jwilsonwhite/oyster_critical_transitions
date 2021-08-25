%will the model reach different steady state conditions when it has
%different starting conditions?

%toggling - the amount of shell and the number of predators.  choose
%whichever way you want to toggle them and they will be improved in the
%future.

% Do EPR calculations using IPM for GTM NERR project
% Started by JW White, April 2019, Edited LS Storch March 2020

%the outer time step in this function is 6 months (because
%that's the timeline of oyster spawning) and the inner time step is 1 week
%(because that's the frequency of environmental data)

%April 2021 - edited input so that it's the initial condition of the oyster
%population, hold all else the same

function [EPR, SSD, P, N, R] = EPR_CriticalTransition(init, addshell, predpenalty, Params, Pred, F, Sal, Temp)%, Pinit)%, Ninit, Sinit, Pinit)%, S) %test where we make the available shell an input
%
% %DEBUGGING
% clear all
% close all
% %DEBUGGING

% %temporary adding christian's salinity data, delete later
% %load('EPR_CT_Params')
% %load('Sal_Christian')
% load('Christian_temp')
% load('Salinity_created_timeseries');
% 
% % %debugging, function part commented out, including function inputs here
% load('EPR_CT_Params')
% % load('testinit')
% % %init = initiwant; %using a "normal" size distribution to start out with
% %init = ones(Params.veclength,1);
% load('testnormalinit')
% init = normaldist;
% addshell = 0;
% pedpenality = 0;
% % %end debugging

%parameters have all been removed from the inside of the function and
%changed to inputs instead.  The file with all of the parameters is
%EPR_CT_Params.m and it produces a .mat file with all of the oyster and
%predator parameters, also named EPR_CT_Params.

%define some more predator things based on input parameters - these are in
%the params file as they should be
% PrefMat = normpdf(repmat(Params.x(:)',[length(Pred.x),1]),...
%     (39.33 + 0.367.*repmat(Pred.x(:),[1,length(Params.x)])),235^0.5);
% % Should sum to one for each predator size class:
% PrefMat = PrefMat./repmat(sum(PrefMat,2),[1,size(PrefMat,2)]);
% Pred.aP = 3.8041.*PrefMat; % attack rate (size dependent)


% %DEBUGGING
% meansal = [ExtremeSummerClimatology, ExtremeSummerClimatology];
% %creating a weird set of sals in which you only go over 30 a few times
% %meansal(meansal < 38) = 20;
% %meansal = meansal*0;
% %meansal(:) = 30; %really bad sal - can we kill everything?
% %meansal(:) = 20;
% Sal = meansal;
% meantemp = [TS_temp, TS_temp, TS_temp, TS_temp];
% %meantemp = meantemp*0;
% %meantemp(:) = 23.2; %hottest average temp col
% Temp = meantemp;
% %END DEBUGGING

%
% Sal = [ExtremeSummerClimatology, ExtremeSummerClimatology]; %TEMPORARILY
% %COMMENTED OUT
% % %Sal = [Sal2(:,1:40), Sal2(:,1:40)];
% % %Sal = Sal2; %the salinity you had in there still wasn't exactly the same as the one in christian's code
% % %Sal = Sal(:,:,1);
% % %Sal = Sal*0+15; %TEMPORARY ADD, DEBUGGING, MAKE SAL CONSTANT
% Temp = [TS_temp, TS_temp, TS_temp, TS_temp]; %New sal timeseries is
% %longer, double time. TEMPORARILY COMMENTED OUT
% % %Temp = TS_temp;
% % %end temporary add

S = nan(1,Params.T); % reef structure
TotBio = nan(1,Params.T); %total biomass, S(t) + living oysters N(t)
N = nan(Params.veclength,Params.T); % oyster population size matrix (semi-annual time steps, based on incoming parameters M and k)
%P = nan(Params.veclength,Params.T); %predator population size matrix
P = nan(Pred.veclength,Params.T); %testing, March 2021
R = nan(Params.veclength,Params.T); % recruitment matrix
Harv = nan(Params.veclength,Params.T); % oysters harvested (semiannual), currently unused

% Initialize the size distributions - only use when debugging and function
% part is commented out!
%N(:,1) = 1; %uncommenting for now, DEBUGGING
%addshell = 0;
%predpenalty = 0;

N(:,1) = init; %initial population distribution is now an input

P(:,1) = 3; %testing, matching with christian's
S(1) = 1;
TotBio(1) = 1;
R(:,1) = 1;

%kernels
kmatNr = kernmatSimp(Params.x,Params,0,'fecundity'); % get fecundity kernel
kmatNg = kernmatSimp(Params.x,Params,0,'growth'); % get growth kernel
kmatNm = kernmatSimp(Params.x,Params,0,'mortality');  % get natural mortality kernel
kmatNf = kernmatSimp(Params.x,Params,F,'fishing');  % get fishing mortality kernel
kmatPr = kernmatSimp(Pred.x,Pred,0,'fecundity'); % predator fecundity kernel
kmatP = kernmatSimp(Pred.x,Pred,0,'both'); %predator growth & mortality
kmatPg = kernmatSimp(Pred.x,Pred,0,'growth'); % predator growth - not currently used just testing March 2021
kmatPm = kernmatSimp(Pred.x,Pred,0,'mortality'); % predator mortality kernel, currently not used just testing March 2021

RunningMean = Sal(1,1); % initialize variable to keep track of running mean of salinity

%Rsalstep = []; %debugging param
%Rshellstep = []; %debugging param
%Rkernelstep = []; %debugging param
Survival = [];
density = [];
A = [];
%Ysave = [];
%NT2SAVE = [];
%NTSAVE = [];

idx = 1;
%howmanysal30s = 0; %want to see how many times the special kernel is activated (total possible iters is Params.T*26

for t = 2:Params.T %outer loop is 6 months because oysters have 2 reproductive seasons in a year
    
    %if t >= 101; keyboard; end
    %June 2021 - calcuating area in a different way
    %density from Byers et al 2015%
    % addshell = addshell/8000; %the ballpark avg density of oysters per 0.25 m^2 is 8kg/0.25 m^2 for southern latitudes (so had to do 8000 g because everything else was in g)
    % addshell = addshell*4; %converting 1/4 m^2 to m^2
    % Units of A need to be m^2, and use the Byers et al conversion.
    
    A(t) = (TotBio(t-1) + addshell(t))/Params.MeanOysterDens*4; % total biomass available for settlement in a timestep (in g) %edited June 2021 to include live oysters in TotBio term.  multiplied by 4 because the mean oyster density is for 0.25m^2
    %A = 0.01; %tester value
    %A = 1000; %tester value
    %added June 2021 - getting an avg area per 0.25 m ^2 based on avg
    A(t) = min(A(t),Params.MaxArea); % limit on total area
    
    %A = (A/Params.density).^2/3; % convert biomass (stored in S) to surface area (in cm^2)
    %A = 100*A; % convert cm^2 to mm^2
    
    %Prey
    Ninter = N(:,t-1); %took out growth in outer loop because it's a weekly rate, only reproduction happens in 6 month loop
    Rinter = kmatNr*Ninter*Params.dx; %intermediate step in R "Rinter" - reproduction
    
    %Rkernelstep(:,t) = Rinter; %debugging
    
    %May 2021 - new salinity based mortality from Lough 1975
    if t == 1 %currently this is never triggered because t starts at 2
        Savg = Sal(1,1); %-M_rs2;
    else
        Savg = nanmean(Sal(23:end,t-1)); %mean of prior 4 weeks
        Tavg = nanmean(Temp(23:end,t-1)); %mean of prior 4 weeks
    end
    
    Sal_save(t) = Savg;
    %percentage surviving over 2 days:
    Y = Params.b0 + Params.b1*Tavg + Params.b2*Savg + Params.b3*Tavg^2 + Params.b4*Savg^2 + Params.b5*(Tavg*Savg);
    Y = Y/100; %it's a % survival but the number is on O(10^1) instead of decimal
    Y(Y<0) = 0; %Y can be negative so have to enforce zero cutoff
    %convert survival to mortality and make daily rate not 2 day rate:
    %M = -log(Y)/2;
    
    %Ysave(t) = Y; %debugging
    %edited june 2021 - applying the salinity based survival and density
    %based survival on the total recruit population as a whole, NOT on
    %individual size classes.  then after salinity and density based
    %survival is applied, re-shape back into a vector so all size classes
    %of recruits are subject to the same mortality
    %Rsum = sum(Rinter);
    %Rsum = Rsum * Y;
    %Y = 0.8368; %DEBUGGING, TEMPORARY
    Rinter = Rinter * Y;%exp(M);
    
    %Rsalstep(:,t) = Rinter; %debugging
    
    %debugging
    %sum(Rinter)
    %pause
    
    %added April 2020 - salinity dependent larval mortality - this totally
    %wipes out the Rinter vector to a comical degree, but all of the
    %parameters were taken from Christian's code, which is Apalachiacola so
    %I don't know what they're supposed to be
    %    %old salinity based mortality
    %     if t == 1 %currently this is never triggered because t starts at 2
    %         Sal_larv = Sal(1,1); %-M_rs2;
    %     else
    %         Sal_larv = nanmean(Sal(23:end,t-1)); %-M_rs2; % mean of prior 4 weeks
    %     end
    %
    %     Sal_save(t) = Sal_larv;
    %
    %     M_r_total = Params.M_r + Params.M_rs*(Sal_larv - Params.M_rs2).^2;
    %     M_r_save(t) = M_r_total;
    
    % end added April 2020
    
    %debugging
    %assignin('base','Rinter',Rinter);
    %sum(Rinter)
    %pause
    
    %old salinity based mortality terms
    %Rinter = Rinter * exp(-M_r_total);
    %Rsave_exp(t) = sum(Rinter);
    
    %sum(Rinter)
    %pause%*Params.dx);
    %remove line immediately below
    %Rinter = Rinter ./ (1 + Rinter.*(pi.*Params.Rmean.^2)./repmat(A,[size(Rinter,1),1])); % survival depends on availability of habitat
    
    %Rshellstep(:,t) = Rinter;
    
    %sum(Rinter)
    %pause
    %R(:,t) = Rinter(:) ./ (1 + Params.DD.*Rsum); % density-dependent mortality, based on Puckett & Eggleston
    %R(:,t) = Rinter ./ (1 + Params.DD.*Rinter); % density-dependent mortality, based on Puckett & Eggleston - %change to later Puckett paper - survival depends linearly based on # of spat
    %added June 2021 - changed density dependence to linear equation from
    %Puckett & Eggleston 2012
    %R(:,t) = Rinter.*(0.7 - 10^(-4).*(Rinter/A));
    Rsum = sum(Rinter);
    density(t) = (Rsum/A(t));
    % Old way: negative exponential. Led to Ricker-like dynamics
    % Survival(t) = exp(-0.00051 * (Rsum/A(t))); %what is the survival rate for the current crop of recruits, based on total density?
    % Survival(t) = max(Survival(t),1e-3);
    % Now: Bev-Holt
    Survival(t) = Params.BH1/(1+Params.BH1/Params.BH2*density(t)); %TEMP COMMENTED OUT DEBUGGING
    Rinter = Rinter.*Survival(t); %convert back to a vector, hitting every size class with the same survival rate
    %Rinter = Rinter.*(0.7 - 10^(-4).*(Rsum/A)); %removed linear
    %relationship because it mapped too much of the survival to zero,
    %replaced with negative exponential above
    % Rinter(Rinter<0) = 0; %have to reinforce zero bound since values are mapped to negative numbers
    R(:,t) = Rinter;
    N(:,t) = Ninter + R(:,t); %final update to population distribution in outer loop.  outer loop only has growth and reproduction
    %Rshellstep(:,t) = Rinter; %debugging
    
    %Predator stuff, stolen from Christian's code
    RR_P = kmatPr*P(:,t-1)*Pred.dx;
    LEP = Pred.LEP; % drill lifetime egg production (R0)
    a = 1/(0.1*LEP); % a = slope at the origin of the Beverton Holt curve
    b = 15; % max density drills = 15/m^2 (from Kimbro 2017 paper)
    RR_P = a*RR_P./(1+(a*RR_P)/b); % pred. DD, Beverton Holt like equation %THE BACKSLASH MAKES THIS A MATRIX IF YOU DON'T HAVE DOTS - WHY DID HE NOT HAVE DOTS?  I PUT IN DOTS, MIGHT CAUSE PROBLEMS LATER, NOT SURE
    P(:,t) = P(:,t-1) + RR_P(:);
    
    %debugging
    %     assignin('base','Rinter',Rinter);
    %     assignin('base','R',R);
    %     assignin('base','N',N);
    %     assignin('base','P',P);
    %     assignin('base','RR_P',RR_P);
    
    
    % Define temporary variables for inner weekly loop
    Nt = N(:,t);
    Pt = P(:,t);
    St = S(t-1);
    % TotBiot = TotBio(t-1); % Not necessary since in the inside loop
    % TotBiot is recalculated each step based on the current values of St
    % and Nt
    
    %Inner weekly loop with temp and salinity information.  Modified from Christian's code
    
    for tt = 1:26 % sub-loop on a weekly timescale for growth & predation dynamics
        
        % Get salinity-dependent kernels if needed - if salinity is really low or really high, increase mortality and inhibit growth
        if Sal(tt,t) <= 5 || Sal(tt,t) >=30
            
            %howmanysal30s = howmanysal30s + 1;
            Params2 = Params; % new params = old params
            Params2.k = Params2.k*0.1; % but no growth
            %Params2.k = 0;
            Params2.M = Params2.M*2; % higher mortality M (100% increase) (arbitrary increase)
            Params2.Mjuv = Params2.Mjuv*2; % higher juv. mort. Mj (100% increase) (arbitrary increase)
            kmatNg = kernmatSimp(Params2.x,Params2,0,'growth'); % use growth kernel with new params
            kmatNm = kernmatSimp(Params2.x,Params2,0,'mortality'); % use natural mortality kernel with new params
            %SAVENG = kmatNg; %debugging
            %SAVENM = kmatNm; %debugging
        else
            kmatNg = kernmatSimp(Params.x,Params,0,'growth'); % use growth kernel with new params
            kmatNm = kernmatSimp(Params.x,Params,0,'mortality'); % use natural mortality kernel with new params
        end
        
        % drill salinity-dependent condition based on Butler 1985 (pg. 5)
        % 10ppt reasonable/consensus cutoff for higher mortality
        if Sal(tt,t) <= 10
            Pred2 = Pred; % new params = old params
            Pred2.M = Pred2.M*2; %0.005; % but higher mortality M (100% increase)
            Pred2.Mjuv = Pred2.Mjuv*2; %0.005; % and higher juv. mort. Mj (100% increase)
            kmatPm = kernmatSimp(Pred2.x,Pred2,0,'mortality'); % use natural mortality kernel with new params
        else
            kmatPm = kernmatSimp(Pred.x,Pred,0,'mortality');  % use regular natural mortality kernel
        end
        
        
        %Now incorporate the various types of mortality and weekly growth info
        Nt2 = kmatNm*Nt; %natural mortality
        
        %NT2SAVE(:,idx) = Nt2; %DEBUGGING
        
        DeadM = Nt - Nt2;
        %Ntf = WkmatNf*Nt2; % fishing, none currently
        %DeadF = Nt - Ntf; %DeadF not used anywhere, no fishing
        Nt = kmatNg*Nt2; %growth - IS THIS SUPPOSED TO BE NT2 OR NT???
        
        %NTSAVE(:,idx) = Nt; %DEBUGGING
        idx = idx+1;
        
        Nt = Nt*Params.dx; %"integrate" at the end when done accounting for growth and mortality (no reproduction in inner loop)
        
        % Predators may change behavior or die due to salinity
        %RunningMean = Sal(1,1); % initialize variable to keep track of running mean of salinity
        
        RunningMean = nanmean([RunningMean,Sal(tt,t)]); % update running mean salinity
        salP = predator_salinity_penalty(Sal(tt:end,t),Temp(tt,t),RunningMean);
        
        % this doesn't seem to be calculating a running mean...
        % it appears to take average of current mean and current salinity
        % not average up to that point
        
        % Predation depends on salinity, habitat complexity, prey density
        aP = Pred.aP ; % attack rate
        hP = Pred.hP ; % handling time % hP is a nxp matrix giving the handling time for pred size p on prey size
        cP = Pred.cP; % Predator interference:
        Tp = Params.Timestep.*salP; % penalize total predation time
        
        % Crawley-Martin function:
        % Crowley-Martin functional response is best fit to mesocosm data.
        % Equation: dN/dt = a*N*P/(1 + b*N + c*P + b*c*N*P)
        %Params.dx*ones(1,length(Params.x)) is replacing Meta.IPM.Pred.Sy because we are no
        %longer using Simpson's because it leaves weird artifacts
        AttackRate = (( Tp*aP'*Pt)) ./ (1 + hP.*Nt + cP.*repmat(Pred.dx*ones(1,length(Pred.x))*Pt,[size(Nt,1),1]) + hP.*cP.*repmat(Pred.dx*ones(1,length(Pred.x))*Pt,[size(Nt,1),1]).*Nt);
        DeadP = AttackRate.*Nt;
        Nt = Nt - DeadP;
        
        % Predators grow & die (depends on salinity)
        %Pt = kmatP*Pt*Pred.dx; %need to do growth and mortality separately
        %if mortality depends on salinity
        Pt = kmatPm*Pt;
        Pt = kmatPg*Pt.*Pred.dx;
        
        % Oysters create new substrate, substrate erodes.
        St = St*exp(-Params.lambdaTAF) + Params.dx*(Params.LW_afdw(:)'*(DeadM + DeadP));
        
        %Added June 2021 - estimate TOTAL biomass - deadshell + living
        %oysters, in grams:
        TotBiot = St + sum(Nt.*Params.LW_afdw(:)); %using AFDW isn't strictly accurate since live oysters are usually measured "wet" but this will just be off by some constant
        
        
        %add updated estimate of total area - St plus living biomass area
        %(current size distribution vector times AFDW) - average biomass
        %from Byers was 8 / 0.25 m^2
        
        %if t == 101 && tt > 1;  keyboard; end
    end % end inner weekly loop over tt
    
    %store results in the semi-annual data
    TotBio(t) = TotBiot;
    N(:,t) = Nt;
    P(:,t) = Pt;
    S(t) = St;
    AR(:,t) = AttackRate;
    
end %end outer 6 month loop

% Stable size distribution
SSD = N(:,end);

% Calcluate EPR (relative, at this point) (this is just biomass at this
% point, fix to be eggs)
EPR = sum(Params.Fec(:).*SSD(:));

%DEBUGGING
% figure(1)
% subplot(1,3,1)
% pcolor(Rkernelstep(1:15,:)) %zooming in on smaller sizes
% ylabel('Size')
% xlabel('Time')
% shading flat
% colorbar
% title('R after kernel step')
% set(gca,'fontsize',16)
% subplot(1,3,2)
% pcolor(Rsalstep(1:15,:)) %zooming in on smaller sizes
% ylabel('Size')
% xlabel('Time')
% shading flat
% colorbar
% title('R after sal step')
% set(gca,'fontsize',16)
% subplot(1,3,3)
% pcolor(Rshellstep(1:15,:)) %zooming in on smaller sizes
% ylabel('Size')
% xlabel('Time')
% shading flat
% colorbar
% title('R after shell dep')
% set(gca,'fontsize',16)
% set(gcf,'color','white')
% 
% figure(2)
% pcolor(N(:,:)) %zooming in on smaller sizes
% ylabel('Size')
% xlabel('Time')
% shading flat
% colorbar
% title('N')
% set(gca,'fontsize',16)
% set(gcf,'color','white')
% 
% figure(3)
% pcolor(R(:,:)) %zooming in on smaller sizes
% ylabel('Size')
% xlabel('Time')
% shading flat
% colorbar
% title('Rshellstep')
% set(gca,'fontsize',16)
% set(gcf,'color','white')
%

%percentspecial = howmanysal30s/(Params.T*26)


end %TEMPORARILY COMMENTED OUT FUNCTION PART