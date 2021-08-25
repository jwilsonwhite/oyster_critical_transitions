%August 2021 - WARNING - the growth kernel is very poorly resolved for
%growth rates near zero, as well as k = 0.  if using a k near zero the
%kernel will create exponential increase of the population instead of
%supressing growth

%WARNING - there are a LOT of hardcoded variables in this code taken from
%oyster_PP_params

%edited january 2020 because now there are separate M and Mjuv so just
%changed the Mj = Param.M line to Mj = Param.Mjuv 

function kxy = mkkernLS(x,y,Param,M_master,Type)

% Construct kernels for IPM
% NOTE: at present there is no salinity effect on growth/mortality, so
% input argument 'Sal' goes unused here

% extract parameters
isjuv_mean = 70; %where do these numbers come from?
isjuv_sd = 2; %where do these numbers come from?

%DEFINE WHO CAN GET FISHED
%isjuv = 1 - normcdf(x,fixparm(4),diff(x(1,1:2))/2);
isjuv = 1 - normcdf(x,isjuv_mean,isjuv_sd); %not currently used

%isjuv = 1 - normcdf(x,Param.Lf,Param.Lfs); % which sizes get fished
isjuv2 = 1 - normcdf(x,Param.j_length,1);  % which sizes experience juvnile mortality rate

%SURVIVAL PART OF KERNEL
switch Type
    case {'mortality','growth','both','fecundity'}
        % natural mortality only
        Ma = Param.M; %adult mortality
        Mj = Param.Mjuv; %juvenile mortality
        m = ones(size(x)).*Ma.*(1-isjuv2) + ones(size(x)).*Mj.*(isjuv2);
        %M = Param.M;
        %m = ones(size(x)).*M; + (1-isjuv).*F; %this is a matrix size x,
        %mortality for each size
    case {'fishing'}
        M = M_master;
        m = (1-isjuv).*M;
end
%this is a matrix size x,
%mortality for each size

p1 = exp(-m); % convert mortality rate to survivorship

%GROWTH PART OF KERNEL
vBs = 0.25; % CV of spread ( std of residuals/length) %taken from oyster_PP_params
growthvar = 0.1; %????

Linf = Param.Linf;
k = Param.k;
%growth
pmean1=Linf - (Linf - x).*exp(-k); % (do not add in x0 for the one-step growth)
%add variability around von Bertalanffy growth
psig1 = pmean1*vBs; % make this a parameter; %.*varparm(3);
%psig1 = pmean1*growthvar; %EPR_NERR has growthvar, oyster_IPM_SS has vBs.. which one???
%evaluate growth part of kernel
p2 = normpdf(y, pmean1, psig1);


%define pr(reproductive) using a maturity ogive
ismat = normcdf(x,Param.Mat,diff(x(1,1:2))/2);
%size(ismat)

pmean2 = repmat(Param.Fec(:)',[length(x),1]);
%size(pmean2)

pmean2 = ismat.*pmean2;
pmean3 = repmat(Param.Rvec(:),[1,length(x)]);
%size(pmean3)

p3 = pmean2 .* pmean3;

p1 = max(0,p1);    %to make sure no negatives
p2 = max(0,p2);    %to make sure no negatives

%added august 2021, the sum of cols don't add up to 1/dx so enforcing that
for i = 1:length(p2)
    p2(:,i) = (p2(:,i)/sum(p2(:,i)))*(1/diff(x(1,1:2)));
end


if k == 0 %trying this out, when k = 0 the kernel goes "funny" so making a special case zero growth kernel where the chances of staying in your current size class are 100%
    p2 = eye(size(p2))*(1/diff(x(1,1:2)));
end
    

switch Type
    case 'growth'
        kxy = p2;
    case {'mortality','fishing'}
        kxy = diag(diag(p1));
    case 'both'
        kxy = p1.*p2;
    case 'fecundity'
        kxy = p3;
end


