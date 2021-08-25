%Aug 2021 - edited such that addshell is now a vector of length t so we can
%decide to only add shell at certain discrete times 

%August 2021 - new section.  basically "done" debugging for now.  we know
%that if we run the code with bad salinity for an extended period of time
%it kills the population.  now if we start with the nearly gone population
%can we rebuild if we have "bad" summers (salinity 30) but normal winters?
%do we have to manually add shell?

clear all
close all

%just creating a bad time series of sal = 30, recreated from inside
%EPR_Criticaltransitions function when i was debugging (now commented out).
% temp is just max hot temp 
load('Christian_temp')
load('Salinity_created_timeseries');
load('EPR_CT_Params')
load('testinit')
load('testnormalinit')

meansal = [ExtremeSummerClimatology, ExtremeSummerClimatology];
meansal = meansal*0;
meansal(:) = 30; %really bad sal - can we kill everything?
%meansal(:) = 20;
Sal = meansal;
meantemp = [TS_temp, TS_temp, TS_temp, TS_temp];
meantemp = meantemp*0;
meantemp(:) = 23.2; %hottest average temp col
Temp = meantemp;
init = ones(Params.veclength,1);
%init = normaldist;
F = 0;
addshell = zeros(1,Params.T);
%what if we only do a few addshells?  
%addshell(round(Params.T/2)) = 100000;
%addshell(round(Params.T/4)) = 100000;

%what if you start adding shell in the middle
for i = round(Params.T/2):Params.T
    addshell(i) = 1000;
end

[EPRold, SSDold, Pold, Nold, Rold] = EPR_CriticalTransition(init, addshell, 0, Params, Pred, F, Sal, Temp);

pcolor(Nold)
shading flat 
colorbar
sum(sum(Nold(:,end)))

%% now what if we try to recover with "normal" winters but still bad summers? (need to run top section first)
load('Christian_temp')
load('Salinity_created_timeseries');
load('EPR_CT_Params')

init = ones(Params.veclength,1);
%newinit = Nold(:,end); %10^-10 population, basically just a negligable recruitment pulse
temp = [TS_temp, TS_temp, TS_temp, TS_temp]; %just putting in a regular temp timeseries for now 
sal = [ExtremeSummerClimatology, ExtremeSummerClimatology]; %based on temp looks like it's summer winter summer winter...
addshell = zeros(1,Params.T);
%how much shell has to be added to significantly affect recovery time?
for i = 1:10:Params.T
    addshell(i) = 100000;
end
    
predpenalty = 0;

for i = 1:2:length(meansal) %setting summers only to high sal 
    sal(:,i) = 30;
end
[EPR, SSD, P, N, R] = EPR_CriticalTransition(init, addshell, predpenalty, Params, Pred, F, sal, temp); %bad summers and bad initial population 

pcolor(N)
shading flat 
sum(sum(N(:,end)))
colorbar

%%  what if you have half of the sal timeseries as bad and half as normal 
load('Christian_temp')
load('Salinity_created_timeseries');
load('EPR_CT_Params')

meansal = [ExtremeSummerClimatology];
meansal = meansal*0;
meansal(:) = 30;
sal = [meansal, ExtremeSummerClimatology]; %half the time series is bad half is normal


init = ones(Params.veclength,1);
%newinit = Nold(:,end); %10^-10 population, basically just a negligable recruitment pulse
temp = [TS_temp, TS_temp, TS_temp, TS_temp]; %just putting in a regular temp timeseries for now 
%sal = [ExtremeSummerClimatology, ExtremeSummerClimatology]; %based on temp looks like it's summer winter summer winter...
addshell = zeros(1,Params.T);
%how much shell has to be added to significantly affect recovery time?
for i = 1:10:Params.T
    addshell(i) = 100000;
end
    
predpenalty = 0;

% for i = 1:2:length(sal) %setting summers only to high sal 
%     sal(:,i) = 30;
% end
[EPR, SSD, P, N, R] = EPR_CriticalTransition(init, addshell, predpenalty, Params, Pred, F, sal, temp); %bad summers and bad initial population 

pcolor(N)
shading flat 
sum(sum(N(:,end)))
colorbar
%% April 2021

%step 1 - if you put in different initial conditions, do you get different
%answers out?  
clear all
close all

load('EPR_CT_Params')

init = ones(Params.veclength,1); %original init in the code is just a vector of ones 

[EPR, SSD, P, N] = EPR_CriticalTransition(init, 0, 0, Params, Pred, F); %'baseline' for the population.  using the SSD as the init now but at different orders of mag

%ok now use the SSD at various orders of mag to see if the model produces
%different outputs based on input 

xtimes = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000];
SSDcell = cell(1,length(xtimes));
EPRcell = cell(1,length(xtimes));
Pcell = cell(1,length(xtimes));
Ncell = cell(1,length(xtimes));
Rcell = cell(1,length(xtimes));

for i = 1:length(SSDcell)
    
    [EPRout, SSDout, Pout, Nout, Rout] = EPR_CriticalTransition(init*xtimes(i), 0, 0, Params, Pred, F);
    EPRcell{i} = EPRout;
    SSDcell{i} = SSDout;
    Pcell{i} = Pout;
    Ncell{i} = Nout;
    Rcell{i} = Rout;
    
end
    
figure(1)
for i = 1:length(SSDcell)
    subplot(5,3,i)
    plot(Params.x,SSDcell{i},'linewidth',2)
    set(gca,'fontsize',16)
    xlabel('size')
    ylabel('abundance')
    title(strcat('init = SSD*',num2str(xtimes(i))));
end
set(gcf,'color','white')

%% if we add a lot of shell does it affect the population?

load('EPR_CT_Params')

init = ones(Params.veclength,1); %original init in the code is just a vector of ones

addshell = [1E-2 1E-1 1E1 1E2, 1E3, 1E4, 1E5, 1E6, 1E7, 1E8, 1E9];

SSDcellshell = cell(1,length(addshell));
EPRcellshell = cell(1,length(addshell));
Pcellshell = cell(1,length(addshell));
Ncellshell = cell(1,length(addshell));
Rcellshell = cell(1,length(addshell));

for i = 1:length(addshell)
    
    [EPR, SSD, P, N, R] = EPR_CriticalTransition(init, addshell(i), 0); %'baseline' for the population.  using the SSD as the init now but at different orders of mag
    
    SSDcellshell{i} = SSD;
    EPRcellshell{i} = EPR;
    Pcellshell{i} = P;
    Ncellshell{i} = N;
    Rcellshell{i} = R;
    
end

figure(2)
for i = 1:length(SSDcellshell)
    subplot(3,4,i)
    plot(Params.x,SSDcellshell{i},'linewidth',2)
    set(gca,'fontsize',16)
    xlabel('size')
    ylabel('abundance')
    title(strcat('added shell = ',num2str(addshell(i))));
end
set(gcf,'color','white')








% %March 2020
% 
% %We want to test the effects of toggling shell and the number of predators
% 
% %grabbing the SSD and shell from the end of a run, this requires running the model to ramp up all of the variables.
% clear all
% close all
% %load('testoutput.mat')
% load('EPR_CT_Params.mat')
% load('EPR_CT_Sal_and_Temp.mat')
% % Sinit = S(end);
% % Ninit = SSD;
% % Pinit = P(end);
% 
% 
% %given these stable values, how many predators can we add before the
% %population crashes?
% 
% iters = 20; %how many times do we want to double the predator pop?
% Ototpop = zeros(1,iters); %save all of the outputs for oysters
% EPRvec = zeros(1,iters);
% SSDmat = cell(1,iters);
% Pcell = cell(1,iters); %save all the predator distributions
% Ptotpop = zeros(1,iters); %save total predator pop
% F = 0; %fishing
% 
% %for i = 1:iters
%   
% %if i > 1
% %   Pinit = Pinit + Pinit; %double the predator each time after the first round
% %end
% 
% %Pcell{i} = Pinit;
% % Ptotpop(i) = sum(Pinit);
% 
% [EPR, SSD, P] = EPR_CriticalTransition(Params, Pred, F, Sal, Temp);%, Ninit, Sinit, Pinit);
% %EPRvec(i) = EPR;
% %SSDmat{i} = SSD;
% %Ototpop(i) = sum(SSD);
% 
% %end
% 
% plot(SSD)
% 
% % plot(Ptotpop,Ototpop,'k')
% % ylabel('Total oyster population')
% % hold on
% % yyaxis right 
% % plot(Ptotpop,EPRvec,'b')
% % legend('Total population','EPR')
% % xlabel('Total predator population')
% % ylabel('EPR')
% % set(gcf,'color','white')
% % set(gca,'fontsize',20)
% % hold off 
% 
% 
% 
% % %script to run EPR_CriticalTransitions with a range of different mortality
% % %values to see the range of solutions
% % 
% % clear all
% % close all
% % load('testingparams.mat')
% % 
% % Params.M = Params.Mjuv; %just putting this in temporarily because when M and Mjuv aren't the same there's a weird dip in the SSD 
% % 
% % %doing a range of 50% above and 50% below the current mortality value just
% % %to see the range of outputs we can get 
% % Mlow = Params.Mjuv*0.5;
% % Mhigh = Params.Mjuv*1.5;
% % 
% % Mrange = linspace(Mlow,Mhigh,100);
% % 
% % EPRvec = zeros(1,length(Mrange));
% % SSDcell = cell(1,length(Mrange));
% % SSDsum = zeros(1,length(Mrange));
% % 
% % for i = 1:length(Mrange)
% %     Params.M = Mrange(i);
% %     Params.Mjuv = Mrange(i);
% %     [EPRvec(i),SSDcell{i}] = EPR_CriticalTransition(Params);
% %     SSDsum(i) = sum(SSDcell{i});
% % end
% % 
% % 
% % figure(1)
% % set(gcf,'position',[100 100 1200 500])
% % subplot(1,2,1)
% % plot(Mrange,EPRvec,'*')
% % xlabel('Mortality rate')
% % ylabel('EPR')
% % title('Mortality rate vs EPR')
% % set(gca,'fontsize',16)
% % subplot(1,2,2)
% % plot(Mrange,SSDsum,'o')
% % xlabel('Mortality rate')
% % ylabel('SSD total')
% % title('Mortality rate vs sum(SSD)')
% % set(gcf,'color','white')
% % set(gca,'fontsize',16)

