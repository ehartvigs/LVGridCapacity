%% ----------------------------------------------------------------------
% This is the main model file for the Swedish version of the network model.
% It contains data loading, processing and saving. None of the data files
% are included in the model package but can be found either on Statistics Sweden or
% the associated publications and projects.
% A description of the model can be found in E. Hartvigsson, et al., "Generating low-voltage grid 
% proxies in order to estimate grid capacity for residential end-use 
% technologies: The case of reisdential solar PV", MethodsX.
% Data can be found in E. Hartvigsson, et al., "Dataset for generating synthetic residential low-voltage
% grids in Sweden, Germany and the UK", Data in Brief
% The model represents parameters and regulations as of 2020.
%% ----------------------------------------------------------------------


clear all
clc
close all


% GIS shapefile can be obtained from Statistics Sweden
filename = 'Inspire_TotPop_Sweref_regionAptShare.shp';
GISData = shaperead(filename,'Selector',{@(v1) (v1<(500000)) && (v1>(1) ), 'Pop'});


% Load annual load profiles from FESTIS project, Swedish Energy Agency, NOT INCLUDED
Demand = load('ResidentialLP.mat');
DemandHH = Demand.HHLoadProfile;
DemandAP = Demand.APTLoadProfile;


% Load solar production dataset, from Z. Norwood, E. Nyholm, T. Otanicar, 
% and F. Johnsson, “A Geospatial Comparison of Distributed Solar Heat and 
% Power in Europe and the US,” PLoS One,  NOT INCLUDED
kk=load('PV_factor_XYcoordinatesSWE.mat');
PV_factor = kk.PV_factor;
PV_X = kk.X;
PV_Y = kk.Y;


% Extract data from the GIS shapefile
x=extractfield(GISData,'X');
x=vec2mat(x,6);
y=extractfield(GISData,'Y');
y=vec2mat(y,6);
[g,tmp] = size(GISData);



pp = 3;     % Voltage, refers to base case
mm = 3;     % Transformer margin, refers to base case
nn = 2;     % Load, refers to base case

%%%%% Control parameters
voltageLimitMatrix = [1.1 0.9; 1.06 0.92; 1.05 0.95];
alphavec = [1.2 1.5 1.8 2.1 2.4];
LoadProfileHH = Load(:,nn);
LoadProfileAP = LoadAP(:,nn);
voltageLimit = voltageLimitMatrix(pp,:);
alpha = alphavec(mm);

% Define variables to save data to

CapacityPerCell = zeros(1,g);
MaxFeeder = zeros(1,g);
Transformers = zeros(1,g);
CapacityPerCustomer = zeros(1,g);
TransformersCap = zeros(1,g);
CustomersPerFeeder = zeros(1,g);
CustomersPerTr = zeros(1,g);
ShareAPT = zeros(1,g);
Limit = zeros(1,g);
PeakLoad = zeros(1,g);
Energy = zeros(1,g);





% Iterate through each cell with data
parfor k=1:g

    % Extract data from GIS dataset
    PopDensity = GISData(k).Pop;
    
    % Find solar PV production profiles for the given cell
    x_cord = mean(y(k,1:5));            % Conversion mixup with X and Y coordinates
    y_cord = mean(x(k,1:5));
    [value ii] = min(abs(PV_X(:)-x_cord) + abs(PV_Y(:)-y_cord));
    [row col] = ind2sub(size(PV_X),ii);
    Sol_factor = PV_factor(row,col,:);
    Sol_factor = squeeze(Sol_factor);
    Sol_factor = repelem(Sol_factor,6);     % Interpolate hourly solar PV production to 10 min
    Sol_factor = Sol_factor(1:52560);
    
    
    [Shareapt customersperfeeder MaxLoad CustomersPerTransformer...
        TrCap LV ll Trans CapPerCust EnergyPerKM2 Limiter]...
        = NetworkModelSweden(PopDensity,LoadProfileHH,LoadProfileAP, Sol_factor, ...
        alpha, voltageLimit);
 
    % Save data from each cell
    CapacityPerCell(k) = ll; 
    MaxFeeder(k) = LV;
    Transformers(k) = Trans;
    CapacityPerCustomer(k) = CapPerCust;
    TransformersCap(k) = TrCap;
    CustomersPerFeeder(k) = customersperfeeder;
    CustomersPerTr(k) = CustomersPerTransformer;
    ShareAPT(k) = Shareapt;
    Limit(k) = Limiter;
    PeakLoad(k) = MaxLoad;
    Energy(k) = EnergyPerKM2;
end


% Save model output
CustomersPerArea = CustomersPerTr.*Transformers;
Energy = double(Energy);
ShareApt = double(ShareAPT);
ModelOutput = GISData;

for h = 1:length(CapacityPerCell)
    ModelOutput(h).Cap = CapacityPerCell(h);
    ModelOutput(h).CapPerCust = CapacityPerCustomer(h);
    ModelOutput(h).Energy = Energy(h);
    ModelOutput(h).Feeder = MaxFeeder(h);
    ModelOutput(h).TrCap = TransformersCap(h);
    ModelOutput(h).TrNumbber = Transformers(h);
    ModelOutput(h).CustPerFeeder = CustomersPerFeeder(h);
    ModelOutput(h).CustPerTr = CustomersPerTr(h);
    ModelOutput(h).Demand = PeakLoad(h);
    ModelOutput(h).ShareAPt = ShareApt(h);
end


filename = sprintf('SWE_DataInBrief2');
shapewrite(ModelOutput,filename);



