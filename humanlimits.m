function ModelOutput = humanlimits(PAR,waterflux,farmed,irrigated,growperiod,eiArray,area,synthetic,lightmove,fusion,contech,longevity,activity,bmi)

% MODEL_OUTPUT is a structure containing information on human population, agricultural
% production, energy use, global resource constraints, and local input constraints ...

% as determined by user inputs to the model, described below:

% Spatial parameters
% PAR is a 3d geospatial array of values representing photosynthetically active radiation (MJ per m2 per day) in each grid cell period
% WATERFLUX is a 3d geospatial array of values representing Precipitation - PET for each grid cell period, in mm or kg per m2
% FARMED is a 2d geospatial array designating the fraction of each grid cell cultivated
% IRRIGATED is a 2d geospatial designating the fraction of each grid cell irrigated
% GROWPERIOD is a binary 3d geospatial array indicating whether or not a given day or composite period falls within the plant growing season
%%% Note: FARMED, IRRIGATED, and GROWPERIOD can be set to a scalar 1, to simplify the specification of extreme scenarios
% EIARRAY is a 3d geospatial array denoting the fraction of light intercepted by crop leaves in each grid cell-day
% AREA is a 2d geospatial array the reports the area (m2) of each grid cell

% Technology parameters
% SYNTHETIC is a binary indicator of whether plants are synthetic and capture additional light in the 700-730 and 750-1075 nm part of the spectrum
% LIGHTMOVE is a binary indicator of whether excess PAR can be stored for use at other times or locations
% FUSION is a binary indicator of whether fusion is viable energy source
% CONTECH is a continuous scalar between 0 and 1, quantifying the level of technology relative to today (0) and the theoretical maximum (1)

% Biological parameters
% LONGEVITY is the deterministic human lifetime
% ACTIVITY is a multiplier on BMR to reflect different levels of daily activity
% BMI is the body mass index ("normal" range is 18.5-25)

% set resource limits, i.e. the maximum availability of several biologically-essential elements (in kg)
resourceLimits(1)= 1*10^20; %Total C stocks in the carbon cycle (sedimentary, atmospheric, terrestrial, marine)
resourceLimits(2)= 8.97*10^22; %Total Ca availability
resourceLimits(3)= 2.48*10^23; %Total Cl availability
resourceLimits(4)= 1.256*10^23; %Total K availability
resourceLimits(5)= 7.714*10^23; %Total Mg availability
resourceLimits(6)= 1.35*10^14 + 3.971*10^18; %Total N availability (terrestrial + atmospheric)
resourceLimits(7) = 1.55*10^23; %Total Na availability
resourceLimits(8)= 5.92*10^21;  %Total P availability
resourceLimits(9)= 1.73*10^23; %Total S availability
resourceLimits(10) = 10665*10^15; %Liquid freshwater availability (kg) in aquifers, lakes, rivers, wetlands, soil, biota

% set parameters
techParams = set_params; % nested function to translate user-supplied tech parameters to model parameters
conversion = techParams(1); % conversion efficiency of intercepted PAR to biomass
harvest = 1 - techParams(2); % harvest index for NPP (i.e. non-principle biomass)
spve = techParams(3); %solar photovoltaic efficiency
resourceLimits(8) = techParams(4)*resourceLimits(8); % P accessible
LMA = techParams(5); % leaf mass per unit area
storage = techParams(6); % energy storage efficiency of batteries
eefp = techParams(7); % production energy required per unit of food energy produced
eip = (1-fusion)/spve*eefp; % solar energy required for production per unit of food energy produced
ei = contech*(1 - eiArray).*growperiod +eiArray; % interception efficiency (only positive during growperiod)

% calculate NPP based on climate and technology parameters
% calculate net energy available for photosynthesis (light within 400-700nm range, unless synthetic plants)
availPAR = (1+synthetic)./(1+0.487*(1+synthetic)*conversion*eip*harvest*ei).*PAR;
% term in denominator uniformly reduces PAR availability by fraction of energy that will be required for agricultural production
% equation assumes light that plants fail to intercept is not available to solar PV

% calculate potential NPP, contingent on growing period and water availability
potentialNPP = conversion*ei.*availPAR; % in MJ per m2 per day

% adjust for growing period and water availability
waterNeed = 0.037663*potentialNPP; % potential water need, in kg per m2 per day
% 0.10809168 kg in 6 mol H2O per 2870 kJ stored energy in 1 mol glucose
% 3.7663e-05 kg per kJ --> 3.7663e-02 kg  per MJ
% keep track of surplus water for later use in identifying local limits
waterDiff = waterflux - growperiod.*waterNeed; % raw water surplus
excessWater = bsxfun(@times,farmed.*(1-irrigated).*area,waterDiff.*growperiod);

[waterfluxIrr, leftover] = allocatewater; % determines temporal allocation of water stored (and delivered via irrigation) within each grid cell
actualNPP = growperiod.*potentialNPP.*max(0,min(1,waterflux./waterNeed)); % without irrigation
actualNPPirr = growperiod.*potentialNPP.*max(0,min(1,waterfluxIrr./waterNeed)); % with irrigation

if lightmove % allow storage and spatio-temporal transfer of PAR for enhanced yield
    [addedNPP,addedNPPirr] = lightmovefcn; % nested function to calculate additional yield
else
    addedNPP = 0;
    addedNPPirr = 0;
end

%convert to harvestable kcal
storedEnergyNonIrr = sum(actualNPP,3).*farmed.*(1-irrigated).*area + sum(addedNPP,3); % sum over time on farmed, non-irrigated areas; total MJ
storedEnergyIrr = sum(actualNPPirr,3).*farmed.*irrigated.*area + addedNPPirr; % sum over time on farmed, irrigated areas; total MJ
storedEnergy = storedEnergyNonIrr + storedEnergyIrr;
kcalNPP = 239.006*storedEnergy; % global net primary product, in kcal (1 MJ = 239.006 kcal)
kcalHarvested = harvest*kcalNPP; % kcal of edible NPP
kcalTotal = sum(sum(kcalHarvested));

humanNeeds = gethumanneeds;
% a vector of per capita flow and stock needs in the following order:
% annual energy (calorie) flow, and embodied stock of C, Ca, Cl, K, Mg, N, Na, P, S, and freshwater

plantNeeds = getplantneeds;

resources4humans = [kcalTotal, resourceLimits] - [0, plantNeeds];

% Apply Liebig's law of the minimum to identify carrying capacity and binding constraint
[CC, global_limit] = min(resources4humans./humanNeeds);

localLimits = findlocal;

ModelOutput = struct('CarryingCapacity',CC,'GlobalLimit',global_limit,'LocalLimits',localLimits,'Harvest',kcalHarvested);


    function techParams = set_params
    % Map the user-supplied technology index value to model parameter values
        best(1) = 0.123; current(1) = best(1)/2.5; % conversion rates of incident PAR during growing season to stored energy, net of respiration
        best(2) = 0; current(2) = 0.5; % residual (non-edible, non-principle) bio frac
        best(3) = 0.868; current(3) = 0.445; %solar photovoltaic efficiency across full light spectrum
        best(4) = 1; current(4) = 2.18*10^13/resourceLimits(8); % fraction of P accessible, based on current proven reserves of P (kg) in phosphate rock
        best(5) = 0.018; current(5) = 0.072; % LMA (kg per m2)
        best(6) = 1; current(6) = 0.9; % battery storage efficiency
        best(7) = 0.000156; current(7) = 0.112; % production energy required per unit of food energy produced
        techParams = contech*(best-current)+current;
    end


    function y = gethumanneeds
    % Determines per capita caloric needs and embodied stock requirements
    % for micronutrients, based on user-supplied values for longevity, bmi
    % and activity levels
        age = 1:longevity;
        height_f = 163.3 - (2*(163.3-151.7))./(exp(0.9261*(age-12.07))+exp(0.1216*(age-12.07))); % based on Model 1 of Preece and Baines (1978) %%% must be in cm for weight and bmr equations
        height_m = 176.9714257 - (2*(176.9714257-163.2799648))./(exp(0.9033929*(age-13.9381565))+exp(0.1007305*(age-13.9381565)));
        weight_f = bmi*(height_f/100).^2; % kg
        weight_m = bmi*(height_m/100).^2;
        bmr_f = 365*(9.99*weight_f + 6.25*height_f - 4.92*age - 161); %kcal per year from Mifflin-St Jeor Equation for females
        bmr_m = 365*(9.99*weight_m + 6.25*height_m - 4.92*age + 5); %kcal per year from Mifflin-St Jeor Equation for males
        bmr = 0.5*bmr_f + 0.5*bmr_m;
        mean_weight = mean(0.5*weight_f + 0.5*weight_m);
        req_kcal = activity*mean(bmr); % per capita energy requirement
        req_stock(1) = 0.18*mean_weight; % C embodied stock requirement
        req_stock(2) = 0.015*mean_weight; % Ca embodied stock requirement
        req_stock(3) = 0.0015*mean_weight; % Cl embodied stock requirement
        req_stock(4) = 0.004*mean_weight; % K embodied stock requirement
        req_stock(5)= 0.001*mean_weight; % Mg embodied stock requirement
        req_stock(6) = 0.03*mean_weight; % N embodied stock requirement
        req_stock(7) = 0.0015*mean_weight; % Na embodied stock requirement
        req_stock(8) = 0.01*mean_weight; % P embodied stock requirement
        req_stock(9) = 0.0025*mean_weight; % S embodied stock requirement
        req_stock(10) = 0.9*0.72*mean_weight; % H2O embodied stock requirement (assuming 10% body fat and 72% water content of lean mass)
        y = [req_kcal, req_stock];
    end

    function z = getplantneeds
    % Determines total stock of each nutrient that must be embodied in plants
        plant_principle = sum(sum(LMA*farmed.*area)); %kg dry mass associated with lowest reported LMA (Poorter et al 2009) necessary to intercept PAR
        plant_stock(1) = 0.45*plant_principle; % kg C
        plant_stock(6) = 0.1761849057*plant_stock(1); % kg N based on Redfield ratio
        plant_stock(2) = 0.083*plant_stock(6); % kg Ca based on Knecht and Goransson (2004)
        %plant_stock(3) = 0; % No data for Cl content
        plant_stock(4) = 0.683*plant_stock(6); % kg K based on Knecht and Goransson (2004)
        plant_stock(5) = 0.087*plant_stock(6); % kg Mg based on Knecht and Goransson (2004)
        %plant_stock(7) = 0; % No data for Na content
        plant_stock(8) = 0.024350441*plant_stock(1); % kg P based on Redfield ratio
        plant_stock(9) = 0.05*plant_stock(6); % kg S based on Inglebeek (2006)
        plant_stock(10) = 3*plant_principle; % Assume 75% of total mass is moisture
        z = plant_stock;
    end

    function [waterflux_irr, leftover] = allocatewater
    % Temporally allocates stored water to its most productive days of use
        noGap = waterflux>=0; % filters out days with negative P-PET, since re-allocating H2O is a waste until total influx > PET
        fluxSum = sum(waterflux.*noGap,3);
        waterDemand = noGap.*growperiod.*waterNeed;
        demandSum = sum(waterDemand,3);
        resid = fluxSum - demandSum;
        resid(resid<0) = 0;
        latentDemand = ~noGap.*growperiod.*(waterNeed - waterflux); %adding absolute value of negative waterflux
        gapFrac = (~noGap.*growperiod.*-waterflux)./latentDemand; % extent of shortfall relative to latent demand
        [~, idxs] = sort(gapFrac,3);
        latentDemSort = zeros(size(idxs)); %preallocate
        lDemMet = false(size(idxs)); %preallocate
        for j=1:size(idxs,1)
            for k=1:size(idxs,2)
                latentDemSort(j,k,:) = latentDemand(j,k,idxs(j,k,:)); % latent demand, sorted by relative shortfall
            end
        end
        cumldem = cumsum(latentDemSort,3);
        lDemMetSort = bsxfun(@le,cumldem,resid); % indicates where reallocated water is sufficient to meet latent demand
        for jj=1:size(idxs,1)
            for kk=1:size(idxs,2)
                lDemMet(jj,kk,idxs(jj,kk,:)) = lDemMetSort(jj,kk,:); % reverses the sort
            end
        end
        demandFrac = bsxfun(@rdivide,waterDemand,demandSum);
        waterflux_irr = bsxfun(@times,fluxSum-resid,demandFrac) + lDemMet.*latentDemand; % adjusted waterflux after reallocation via storage and irrigation
        leftover = resid-sum(lDemMet.*latentDemand,3); % water not usefully reallocated
        leftover = leftover.*area.*farmed.*irrigated; % quantify in absolute (not per m2) units, restricted to relevant areas
    end

    function [addedNPP,addedNPPirr] = lightmovefcn
    % Spatially reallocates excess, stored PAR to its most productive location
        hypotheticalNPP = bsxfun(@times,area.*farmed.*(1-irrigated),potentialNPP - actualNPP);
        hypotheticalNPPirr = bsxfun(@times,area.*farmed.*irrigated,potentialNPP - actualNPPirr);
        hypotheticalSum = sum(sum(sum(hypotheticalNPP))); %NPP-equivalent excess PAR
        hypotheticalSumIrr = sum(sum(sum(hypotheticalNPPirr)));
        excessSum = sum(sum(sum(excessWater))); % excess water during growing season in farmed area
        excessSumIrr = sum(sum(leftover)); % excess water in farmed, irrigated area
        excessFrac = 1/excessSum*excessWater;
        excessFrac(~isfinite(excessFrac)) = 0;
        excessRatio = excessWater./(0.037663*hypotheticalSum*excessFrac);
        excessRatio(~isfinite(excessRatio)) = 0; % eliminate Nan or Inf
        excessFracIrr = 1/excessSumIrr*leftover;
        excessRatioIrr = leftover./(0.037663*hypotheticalSumIrr*excessFracIrr);
        excessRatioIrr(~isfinite(excessRatioIrr)) = 0;
        addedNPP = storage*hypotheticalSum*excessFrac.*min(1,excessRatio);
        addedNPPirr = storage*hypotheticalSumIrr*excessFracIrr.*min(1,excessRatioIrr);
    end

    function localLimits = findlocal
    % Determines NPP forgone due to insufficient PAR vs insufficient water,
    % and reports the fraction of forgone NPP due to lack of water
        excessWaterIrr = bsxfun(@times,farmed.*irrigated.*area,(waterfluxIrr - waterNeed).*growperiod);
        insufficientWater = -excessWater;
        insufficientWaterIrr = -excessWaterIrr;
        excessWater(excessWater<0) = 0;
        insufficientWater(insufficientWater<0) = 0;
        insufficientWaterIrr(insufficientWaterIrr<0) = 0;
        forgoneNPPdue2PAR = (1-irrigated).*(1/0.037663*sum(excessWater,3)-sum(addedNPP,3)) + irrigated.*(1/0.037663*leftover-addedNPPirr);
        forgoneNPPdue2water = (1-irrigated).*(1/0.037663*sum(insufficientWater,3)) + irrigated.*(1/0.037663*sum(insufficientWaterIrr,3));
        forgoneNPP = forgoneNPPdue2PAR + forgoneNPPdue2water;
        localLimits = forgoneNPPdue2water./forgoneNPP;
    end
end