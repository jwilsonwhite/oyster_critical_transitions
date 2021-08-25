function SalP = predator_salinity_penalty(Sal,Temp,RunningMean)

% Penalize predation rate due to low-salinity exposure
% Parameterized based on mesocosm trials

% In Sal, col 1 is weekly salinity; col 2 is #days (of 7) < 10 ppt, col 3 is 5-day < 5 ppt exceeded

if Temp < 20 % Use threshold of 15º for temperature limit 
    
    TempP = max(0,(Temp-15)./20); % penalty as drop below 20 degrees, zero at 15 degrees
    
else
    TempP = 1;
end

% Overall salinity-dependence of predation (based on Kimbro field
% experiments):
SalP1 = max(0,min(1,1-(30.6-Sal(1)).*0.05));

% Reduction due to sudden changes in salinity:

if RunningMean - Sal(1) > 10 % 10 ppt reduction in mean salinity
    SalP2 = 0;
else % if not such a big reduction

Reduction = max(0,RunningMean - Sal(1)); % if salinity goes up, no change

SalP2 = max(0,(1 - 0.1*Reduction)); % based on % reduction observed for 10-day exposures to 5 & 10 ppt drops (updated 8/2016)
    
end

SalP = TempP .* SalP1 .* SalP2;
