%Parameters for the EPR_CriticalTransitions function 
%M, Mjuv, k, Linf for the oysters updated with Params from Christian's code
%(because before I was using params from GTMNERR which is the wrong area)

%Oyster PARAMETERS - a lot of these were taken from oyster_PP_params or EPR_NERR or Christian's code
Params.k = 0.0122; %taken from Christian's code - weekly 
Params.Linf = 120.3500; %taken from Christian's code - mms
Params.Mjuv = 0.0186; %taken from Christian's code - Based on 80% survival over 12 weeks in cages from Kimbro experiment A2ii
Params.M = 0.0052; %taken from Christian's code - Based on 94% survival over 12 weeks in cages from Kimbro experiment A1v
Params.T = 160; %using salinity created from OysterIPM_ExtremeDisturbances so have longer time series now 
%Params.T = 62; %how many times do we iterate to get the stable size distribution? changed to 62 because that's how many cols are in the sal and temp data
Params.Timestep = 1;
Params.veclength = 250; %number of patches in the vector
%Params.x = linspace(0,Params.Linf*2,Params.veclength); % space vector
Params.x = linspace(0,150,Params.veclength); % space vector %TEMPORARILY CHANGED TO MATCH CHRISTIAN'S 
Params.dx = diff(Params.x(1:2)); %spacing, used for integrating
Params.LW_afdw = 5.09e-5*Params.x.^2.365; % Ash Free Dry Weight; units = g ; derived from Kimbro field samples in salinity zone 2 (taken from oyster_PP_params)
Params.Mat = 40; %maturity parameter, taken from oyster_PP_params
Params.Fec = 19.86e6 .* Params.LW_afdw.^(1.17); %fecundity parameter, taken from oyster_PP_params
Params.Fec = Params.Fec .* (Params.x>=35); %Impose size at maturity of 35 mm
Params.Rmean = 2.5; % mean spat size in week 1 - taken from Christian's code 
Params.Rstd = 0.5; % std in spat size - taken from Christian's code 
Params.Rvec = normpdf(Params.x,Params.Rmean,Params.Rstd); %taken from EPR_NERR, initializing recruitment vector
Params.DD = 1e-3; % based on Puckett & Eggleston (2012), Fig 7. (taken from oyster_PP_params)
Params.j_length = 15; %juvenille mortality size limit - taken from Christian's 
Params.isjuv = Params.x < 76.2;
F = 0; %fishing rate
Params.lambdaTAF = 0.1/52;  % Based on annual rate of 0.1 from Powell et al. (2012)
Params.density = 0.849; % grams/cm^3. From aqua-calc.com
Params.MeanOysterDens = 8000; %the ballpark avg density of oysters per 0.25 m^2 is 8kg/0.25 m^2 for southern latitudes (so had to do 8000 g because everything else was in g) Byers et al 2015

Params.MaxArea = 20e6; % Cat Point Bar is ~20 km^2 (probably somewhat less actually)

%salinity based larval mortality parameters - taken from Christian's code 
Params.M_r = 7.8; % 0.26 per day * 30 d (Rumrill 1990)
Params.M_rs = 2.934; % fitted param from SSIPM
Params.M_rs2 = 15; % optimum salinity (Davis 1958 Biol Bull)

%new salinity based mortality coefficients from Lough 1975, 2 day survival
Params.b0 = -643.9149; %constant
Params.b1 = 27.7755; %T
Params.b2 = 32.8617; %S
Params.b3 = -0.5195; %T^2
Params.b4 = -0.6234; %S^s
Params.b5 = -0.0971; %T*S

%July 2021 - Will fit the data to a Beverton-Holt model with the following
%params (this was replacing other density dependence mechanisms that were
%mapping too much to zero) 
Params.BH1 = 0.8254;
Params.BH2 = 2202; 


%predator parameters, all stolen from christian's code - are these weekly, monthly... ?????
%Pred.x = linspace(0,100,Params.veclength); %added April 2020, trying the pred space vector with a different scale than the prey space vector
Pred.veclength = 100;
Pred.x = linspace(0,100,Pred.veclength); %Temp add March 2021, testing
Pred.dx = diff(Pred.x(1:2));
Pred.LEP = 2332400;
Pred.k = 0.1; % units = 1 week; based on info in Butler (1985), see notebook p. 10
Pred.Linf = 60;  % units = mm ****
Pred.M = 0.0025; % units = 1 week. Based on 12% annual loss rate (Butler 1985)
Pred.Mjuv = 0.0025;
Pred.vBs = 0.1; % CV spread for growth
Pred.cP = 0.5203;
Pred.hP = 1.1622;
%Pred.aP = 3.8041;
%Pred.Fec = 43350 * Params.x.^(0.782); % maybe 3e5 .* Meta.IPM.Pred.x.^(0.782); ???
Pred.Fec = 43350 * Pred.x.^(0.782); %changed March 2021, testing 
Pred.Rmean = 2; % mean drill size in week 1
Pred.Rstd = 0.5; % std in spat size
%Pred.Rvec = normpdf(Params.x,Pred.Rmean,Pred.Rstd);
Pred.Rvec = normpdf(Pred.x,Pred.Rmean,Pred.Rstd); %changed March 2021, testing 
Pred.Mat = 25; % 2.5 cm, smallest size observed producing viable embryos, (Butler 1985)
Pred.j_length = 21; %this doesn't matter because Mj and M are the same - source: Butler, P. A. (1985). Synoptic review...southern oyster drill
% Adjust attack rate based on size preference of drills vs oysters.
% (Mesocosm, Report Fig. 23A)
PrefMat = normpdf(repmat(Params.x(:)',[length(Pred.x),1]),...
    (39.33 + 0.367.*repmat(Pred.x(:),[1,length(Params.x)])),235^0.5);
% Should sum to one for each predator size class:
PrefMat = PrefMat./repmat(sum(PrefMat,2),[1,size(PrefMat,2)]);
Pred.aP = 3.8041.*PrefMat; % attack rate (size dependent)

save('EPR_CT_Params.mat')

