%% ------------------------------------------------------------------------
% This is the network generator, the inputs are population density, load
% profile for the customers, solar production profiles, and external model
% parameters found in the MainModelSweden file.
%% ------------------------------------------------------------------------



function [frac_AP,CustomersPerFeeder,MaxLoad,CustomersPerTransformer,TransformerType...
    ,LLVMax,solcap,NumberOfTransformers,CapPerCustomer,EnergyPerKM2,Limiter]...
    = NetworkModelSweden(PopDensity,HouseLoadProfile, AptLoadProfile1, SolProfile, Alpha, VoltageLimit)


Design_Z = 0.65;    % Design maximum earth loop impedance
Design_voltage = [1.06 0.92];   % Design voltage


%%  Set model parameters  %%

% Demographic parameters
PeoplePerHouse = 2.7;
PeoplePerApt = 2;

% Nominal voltage
Vn = 400;

% Fuse ratings
AptFuseRating = 10;
HouseFuseRating = 20;

% Electricity consumption parameters
Limiter = 0;
margin = Alpha;               % Transformer sizing margin
AptBuilding_share = 1.25;     % Increase in electricity consumption due to ancillary electricity usage in apartment buildings
                              % https://www.sciencedirect.com/science/article/pii/S0378778811005019#fig0010
AVGEnergyHouse = 18.5*1000;     % Average annual electricity consumption for detached house in Sweden
AVGEnergyApt = 3.71*1000;       % Average annual electricity consumption for apartment in Sweden
AVGEnergyApt = AVGEnergyApt*AptBuilding_share;
solcap = 0;
PowerFactor = 0.9;
PowerFactor = sind(acosd(PowerFactor));

% Velander coefficients
k2_House = (0.025+0.05)/2;
k1_House = 0.0003;
k2_Apt = 0.014;
k1_Apt = 0.000264;


% Power system parameters
TransformerCap = [1250 800 500 315 200 100 50];
TransformerCost = [195272 134751 101565 70501 53509 38446 32140];
Transformer_Impedance = [6.5 10 13 20 32 65 130]/1000;
Z_TransformerR_list = [0.000848114120254821 0.00148841681507050 0.00251028652976876 0.00492770793724613 0.00776114000116266 0.0202385770250776 0.0404771540501553];
Z_TransformerX_list = [0.00763302708229339 0.0119073345205640 0.0125514326488438 0.0197108317489845 0.0310445600046506 0.0607157310752329 0.121431462150466];

CableCapacity = [52 67 80 94 138 178 230];
fuses_tripping_size = [6 10 16 20 25 32 40 50 63 80 100 125 160 200 230;...     % Table 3 in EN60269, for fuses less than 16A from BS3036
    18 32 65 85 110 150 190 250 320 425 580 715 950 1250 1490];               % 1490 obtained through linear interpolation
targetFuseRatings = [6 10 16 20 25 32 40 50 63 80 100 125 160 200 230];
Z_lineR_list = [1.83 1.15 0.76 0.641 0.32 0.206 0.125];
Z_lineX_list = [0.0861 0.0817 0.0783 0.0779 0.0762 0.0745 0.0752];
targetLineImpedance = [4.18 2.63 1.91 1.47 0.746 0.495 0.324];

% Cable sizing parameters
c = 0.85;
Ufn = 230;


% Cable costs per km (SEK)
if PopDensity>1000    
    CostPerKmLineLV = 827000;
    CostPerKmLineMV = 1140746;
elseif PopDensity>200   
    CostPerKmLineMV = 887790;
    CostPerKmLineLV = 540000;
else 
    CostPerKmLineMV = 512337;
    CostPerKmLineLV = 177000;
end


% Set feeder parameters based on demographic area
if PopDensity<=200         
    LLAF_LV = 1;            % Feeder adjustment factor (gamma)                  
elseif PopDensity>200 && PopDensity<=1000
    LLAF_LV = 1.1;    
else
    LLAF_LV = 1.2; 
end


% Calculate the share of single family and multifamily households.
if PopDensity <= 4000
    frac_AP = min([(-2.5797*10^(-8)*PopDensity^2 + 0.000257*PopDensity + 0.3782) 1]);
    frac_AP = max([frac_AP 0]);
    frac_HH = 1-frac_AP;
else
    frac_AP = 1;
    frac_HH = 0;
end


% Average fuse size used in grid design
fuse = AptFuseRating*frac_AP+HouseFuseRating*frac_HH;

% Household composition
PeoplePerHousehold = PeoplePerHouse*frac_HH + PeoplePerApt*frac_AP;


% Special cases of low population density
if PopDensity<4
    NumberOfCustomers = PopDensity;
    frac_AP = 0;
    frac_HH = 1;
else
    NumberOfCustomers = round(PopDensity/PeoplePerHousehold);
end


% Construct average load profile
AVG_LoadProfile = (frac_HH.*HouseLoadProfile+frac_AP.*AptLoadProfile1);

% Calcualte average electricity consumption
AVGEnergy = (frac_AP*AVGEnergyApt + frac_HH*AVGEnergyHouse);

% Set initial values for cost-minimization 
CustomersPerTransformer = NumberOfCustomers;
PeakLoad = margin*(k1_House*CustomersPerTransformer*frac_HH*AVGEnergyHouse+k2_House*sqrt(CustomersPerTransformer*frac_HH*AVGEnergyHouse)) + (k1_Apt*CustomersPerTransformer*frac_AP*AVGEnergyApt+k2_Apt*sqrt(CustomersPerTransformer*frac_AP*AVGEnergyApt));
n0 = ceil(PeakLoad/TransformerCap(1));
ToTCost = zeros((200),length(TransformerCap));
MaxLoad = PeakLoad/margin;

% Finding lest cost option of transformer capacity and number of
% transformers to supply the load.
for gg = n0:round(NumberOfCustomers)
    
    CustomersPerTransformer = round(NumberOfCustomers/gg);
    %coincidenceTR = (k1_House*CustomersPerTransformer*AVGEnergy+k2_House*sqrt(CustomersPerTransformer*AVGEnergy))/(CustomersPerTransformer*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));
    PeakLoadTR = margin*((k1_House*CustomersPerTransformer*frac_HH*AVGEnergyHouse+k2_House*sqrt(CustomersPerTransformer*frac_HH*AVGEnergyHouse)) + (k1_Apt*CustomersPerTransformer*frac_AP*AVGEnergyApt+k2_Apt*sqrt(CustomersPerTransformer*frac_AP*AVGEnergyApt)));
    
    
    indexTR = find((PeakLoadTR./(TransformerCap)<1)~=0,1, 'last');
    TransformerCostArea = TransformerCost(indexTR)*gg;
    
    % Size of area supplied by on transformer is in km^2
    A_TS = 1/gg;
    NubmerOfConnectionsTransformer = round(NumberOfCustomers/gg); % Number of connections per transformer

    %ConnectionPointDensity = NubmerOfConnectionsTransformer/A_TS;
    
    % Calcualte distance between low (d) and medium (dMV)
    % voltage supply points
    d = sqrt(A_TS)/(sqrt(NubmerOfConnectionsTransformer)+1)+0.02;
    dMV = 1/(sqrt(gg)+1);
   
    if mod(round(sqrt(NubmerOfConnectionsTransformer)),2)
        n = NubmerOfConnectionsTransformer-1;
    else
        n = NubmerOfConnectionsTransformer + sqrt(NubmerOfConnectionsTransformer) -2;
    end
    
    if round(sqrt(gg)) == 1
        nMV = 0.5;
    elseif mod(round(sqrt(gg)),2)
        nMV = gg-1;
    else
        nMV = gg + sqrt(gg) -2;
    end
    
    LengthLVPerTransformer = n*d;
    LengthMV = nMV*dMV;
    MV_cost = LengthMV*CostPerKmLineMV;
    LV_cost = gg*LengthLVPerTransformer*CostPerKmLineLV;
    
    LineCost = LV_cost + MV_cost;
    ToTCost((gg-n0+1),indexTR) = TransformerCostArea + LineCost;
end

% Set any potential 0 costs option to NaN
ToTCost(ToTCost==0)=NaN;

% Find minimum cost that supply the area, and the corresponding combination
% of transformer capacity and numbers
MinCost=min(min(ToTCost));
size(ToTCost);
[xt,TrSize]=find(ToTCost==MinCost);
xt=xt(1);
TrSize=TrSize(1);
NumberOfTransformers = n0+xt-1;
TransformerType = TransformerCap(TrSize);
Z_transformer = Transformer_Impedance(TrSize);

% Update area parameters
A_TS = 1/NumberOfTransformers;

CustomersPerTransformer= round(NumberOfCustomers/NumberOfTransformers);
coincidenceTR = (k1_House*CustomersPerTransformer*AVGEnergy+k2_House*sqrt(CustomersPerTransformer*AVGEnergy))/(CustomersPerTransformer*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));


% Handle special cases of the uniform distribution of customers
if CustomersPerTransformer == 1
    LLVMax = 0.1;
    d = LLVMax;
    n = 1;
    nc = 1;
elseif CustomersPerTransformer == 2
    LLVMax =  LLAF_LV*sqrt(A_TS)/4 + 0.02;
    d = LLVMax;
    n = 2;
    nc = 2;
else
    n = round(sqrt(CustomersPerTransformer));
    nc = round(CustomersPerTransformer/4);
    if n == 2
        d = LLAF_LV*sqrt(A_TS)/(n);
        LLVMax = d;
        nc = 2;
    else
        LLVMax = LLAF_LV*sqrt(A_TS)*(n-1)/(n)+0.02;
    end
end




% Calcualte transformer resistance and reactance
Transformer_R = Z_TransformerR_list(TrSize);
Transformer_X = Z_TransformerX_list(TrSize);


% Set branches of feeders based on number of supplied customers
if n > 10
    m = zeros(1,5);
    mc = zeros(1,5);
    m(1)= round((n-1)/5);
    m(2) = m(1);
    m(3) = m(1);
    m(4) = m(1);
    m(5) = n-m(1)-m(2)-m(3)-m(4);
    mc(1)= round((nc-1)/5);
    mc(2) = mc(1);
    mc(3) = mc(1);
    mc(4) = mc(1);
    mc(5) = nc-mc(1)-mc(2)-mc(3)-mc(4);
    p = 5;
elseif n > 8
    m = zeros(1,4);
    mc = zeros(1,4);
    m(1)= round((n-1)/4);
    m(2) = m(1);
    m(3) = m(1);
    m(4) = n-m(1)-m(2)-m(3);
    mc(1)= round((nc-1)/4);
    mc(2) = mc(1);
    mc(3) = mc(1);
    mc(4) = nc-mc(1)-mc(2)-mc(3);
    p = 4;
elseif n>6
    m = zeros(1,3);
    mc = zeros(1,3);
    m(1)= round((n-1)/3);
    m(2) = m(1);
    m(3) = n-m(1)-m(2);
    mc(1)= round((nc-1)/3);
    mc(2) = mc(1);
    mc(3) = nc-mc(1)-mc(2);
    p = 3;
elseif n == 1 || n == 2
    m = zeros(1,1);
    mc = zeros(1,1);
    m(1) = n;
    mc(1) = nc;
    p = 1;
else
    m = zeros(1,2);
    mc = zeros(1,2);
    m(1)= round((n-1)/2);
    m(2) = n-m(1);
    mc(1)= round((nc-1)/2);
    mc(2) = nc-mc(1);
    p = 2;
end

mc = flip(mc);
d_long = zeros(1,p);

% Calculate length of each feeder section, and aggregate for a total feeder
% length.
for uu = 1 :p
    if n == 1
        d_long(uu) = (m(uu))*d*1000;
    else
        d_long(uu) = (m(uu)-1/p)*d*1000;
    end
end
d_long = flip(d_long);


% Set variables used in cable sizing
Cable = zeros(1,p);
ixCable = [];
R = zeros(1,p);
X = zeros(1,p);
ZmaxThick = zeros(1,p);
coincidenceLine = zeros(1,p);
P_demand = zeros(1,p);
mp = zeros(1,p);
z_earth = zeros(1,p+1);
RX_multiplier = ones(1,p);
L_max = zeros(7,p);
I_line_fuse = zeros(1,p);
Z_Line = targetLineImpedance;
CustomersPerFeeder = sum(mc);
Z_loop = [Z_transformer*1000];


% Size feeders based on fuse ratings and tripping criteria
for rt=1:p
    rr = sum(mc(rt:end));
    mp(rt) = rr;
    coincidenceLine(rt) = (k1_House.*rr.*AVGEnergy+k2_House.*sqrt(rr.*AVGEnergy))./(rr.*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy))); 
    P_demand(rt) = max(AVG_LoadProfile)*rr*coincidenceLine(rt);
    I_line = coincidenceLine(rt).*fuse*rr;
    
    % Check and update fuse ratings based om demand
    while I_line>max(targetFuseRatings)
        rr = round(rr/2);
        RX_multiplier(rt) =  RX_multiplier(rt)*0.5;
        coincidenceLine(rt) = (k1_House.*rr.*AVGEnergy+k2_House.*sqrt(rr.*AVGEnergy))./(rr.*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));
        P_demand(rt) = max(AVG_LoadProfile)*rr*coincidenceLine(rt);
        I_line = coincidenceLine(rt).*fuse*rr;
    end
    
    
    I_line_fuse(rt) = targetFuseRatings(find((ceil(I_line)./targetFuseRatings<=1)~=0,1, 'first'));
    index = find((I_line_fuse(rt)./targetFuseRatings<=1)~=0,1, 'first') ;
    Iu = fuses_tripping_size(2,(index));
    ZmaxThick(rt) = 1000*c*Ufn/Iu;
    
    
    
    
    % Set maximum length of feeders based on available capacity and fuse
    % ratings
    L_max(:,rt) = (ZmaxThick(rt)-sum(Z_loop))./targetLineImpedance';
    while sum(L_max(:,rt)>d_long(rt)) == 0
        rr = round(rr/2);
        RX_multiplier(rt) = 0.5;
        coincidenceLine(rt) = (k1_House.*rr.*AVGEnergy+k2_House.*sqrt(rr.*AVGEnergy))./(rr.*(k1_House*AVGEnergy+k2_House*sqrt(AVGEnergy)));
        P_demand(rt) = max(AVG_LoadProfile)*rr*coincidenceLine(rt);
        I_line = coincidenceLine(rt).*fuse*rr;
        I_line_fuse(rt) = targetFuseRatings(find((ceil(I_line)./targetFuseRatings<1)~=0,1, 'first'));
        index = find((I_line_fuse(rt)./targetFuseRatings<=1)~=0,1, 'first') ;
        Iu = fuses_tripping_size(2,(index));
        ZmaxThick(rt) = 1000*c*Ufn/Iu;
        L_max(:,rt) = (ZmaxThick(rt)-sum(Z_loop))./targetLineImpedance';
    end
    
    
    
    L_max(L_max<d_long(rt))=NaN;
    L_max(CableCapacity<I_line_fuse(rt),rt) = NaN;
    [void ixZloop] = min(L_max(:,rt));
    Z_loopNew = targetLineImpedance(ixZloop)*d_long(rt);
    Z_loop = [Z_loop Z_loopNew];
end


L_max(isnan(L_max))=0;
z_earth(1) = Z_transformer*1000;

% Check cable thermal capacity is within limits
for gg=1:p
    ixCable(gg) = find(L_max(:,gg),1);
    Cable(gg) = CableCapacity(ixCable(gg));
    if CableCapacity(ixCable(gg))<I_line_fuse(gg)
        if I_line_fuse(gg) == 250 || I_line_fuse(gg) == 315
            break
        else
            ixCable(gg) = find((I_line_fuse(gg)./CableCapacity<=1)~=0,1, 'first');
            Cable(gg) = CableCapacity(ixCable(gg));
        end
    end
    R(gg) = Z_lineR_list(ixCable(gg))*d_long(gg)*RX_multiplier(gg);
    X(gg) = Z_lineX_list(ixCable(gg))*d_long(gg)*RX_multiplier(gg);
    z_earth(gg+1) = z_earth(gg) + targetLineImpedance(ixCable(gg))*d_long(gg)*RX_multiplier(gg);
    
end


% Update cable capacity, resistance, reactance and erath impedance
id = find(ixCable==max(ixCable),1);
for hh=1:id
    Cable(hh) = CableCapacity(ixCable(id));
    R(hh) = Z_lineR_list(ixCable(id))*d_long(id);
    X(hh) = Z_lineX_list(ixCable(id))*d_long(id);
    z_earth(hh+1) = z_earth(hh) + targetLineImpedance(ixCable(hh))*d_long(hh)*RX_multiplier(hh);
end




% Convert reistance and reatice to correct quantitiy.
R = R./1000;
X = X./1000;



% Check voltage drop along feeder is within regulations and update if
% voltage is outside limits
V = zeros(1,p+1);
V(1) = Vn - coincidenceTR*CustomersPerTransformer*(Transformer_R.*(fuse*230*3)+Transformer_X.*(fuse*230*3)*PowerFactor)./(Vn);

for kk=1:p
    V(kk+1) = V(kk)-(R(kk).*(P_demand(kk))+(X(kk)).*P_demand(kk)*PowerFactor)./(Vn);
end

while  min(V./Vn)<Design_voltage(2)
    [val pos] = min(ixCable);
    if min(ixCable) == 7
        [val2 pos] = max(RX_multiplier);
        RX_multiplier(pos) = RX_multiplier(pos)*0.5;
        R(pos) = Z_lineR_list(val)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val;
        Cable(pos) = CableCapacity(val);
        z_earth(pos+1) = z_earth(pos) + targetLineImpedance((val))*d_long(pos)*RX_multiplier(pos);
    else
        R(pos) = Z_lineR_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val + 1;
        Cable(pos) = CableCapacity(val+1);
        z_earth(pos+1) = z_earth(pos) + targetLineImpedance((val+1))*d_long(pos)*RX_multiplier(pos);
    end
    
    
    V = zeros(1,p+1);
    V(1) = Vn -min((coincidenceTR*CustomersPerTransformer*(Transformer_R.*(fuse*230*3)+Transformer_X.*(fuse*230*3).*PowerFactor)./Vn));
    for kk=1:p
        V(kk+1) = V(kk)-(R(kk).*(P_demand(kk).*1000)+(X(kk)).*P_demand(kk).*1000*PowerFactor)./(Vn);
    end
    
end



% Checking that loop impedance is lower than maximum value, otherwise
% update cable parameters
z_earth = z_earth/1000;
while  max(z_earth)>Design_Z
    
    [val pos] = min(ixCable);
    if min(ixCable) == 7
        [val2 pos] = max(RX_multiplier);
        RX_multiplier(pos) = RX_multiplier(pos)*0.5;
        
        for gg = pos:length(ixCable)
            z_earth(gg+1) = z_earth(gg) + Z_Line(val+1)*d_long(gg)*RX_multiplier(gg)/1000;
        end
        R(pos) = Z_lineR_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val;
        Cable(pos) = CableCapacity(val)*(1/RX_multiplier(pos));
    else
        for gg = pos:length(ixCable)
            z_earth(gg+1) = z_earth(gg) + Z_Line(val+1)*d_long(gg)*RX_multiplier(gg)/1000;
        end
        R(pos) = Z_lineR_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val+1;
        Cable(pos) = CableCapacity(val+1);
    end
    
    
end

%z_loop = repelem(z_earth,NumberOfTransformers);
%CableSize = mean(Cable);




% Add solar PV systems to each connection point on the feeder in increments
% of 0.5 kW
for ll=0.5:0.5:100
    
    % Calculate new voltage profile
    V0 = 400 - min((CustomersPerTransformer*(Transformer_R.*(coincidenceTR*fuse*230*3-(ll).*SolProfile.*1000)+Transformer_X.*(coincidenceTR*fuse*230*3))./Vn));
    if length(mp) == 1
        voltage = V0 - min((mp.*(R.*(coincidenceLine.*AVG_LoadProfile-(ll).*SolProfile.*1000)+X.*(coincidenceLine.*AVG_LoadProfile*PowerFactor))./Vn)');
    else
        voltage = V0 - min(sum((mp.*(R.*(coincidenceLine.*AVG_LoadProfile-(ll).*SolProfile.*1000)+X.*(coincidenceLine.*AVG_LoadProfile*PowerFactor))./Vn)'));
    end
    
    % Calcualte maximum power level
    deltaCurrentPerCustomerDouble = mc.*max(abs(coincidenceLine.*AVG_LoadProfile-ll.*SolProfile.*1000))/400/sqrt(3);
    deltaPowerPerCustomer = max(abs(coincidenceTR*AVG_LoadProfile-ll.*SolProfile.*1000));
    deltaPower = deltaPowerPerCustomer.*CustomersPerTransformer;
    
    
    % Check against voltage and thermal limits and abort if limits are
    % reached.
    if (max(max(voltage))/400)>VoltageLimit(1)
        solcap = (ll-0.5)*NumberOfCustomers;
        Limiter = 10;
        break
        
    elseif deltaCurrentPerCustomerDouble>Cable
        solcap = (ll-0.5)*NumberOfCustomers;
        Limiter = 25;
        break
        
    elseif deltaPower>(TransformerType*1000)
        solcap = (ll-0.5)*NumberOfCustomers;
        Limiter = 20;
        break
        
    end
end


EnergyPerKM2 = (NumberOfCustomers)*(ll-0.5)*sum(SolProfile)/6;       % i kW
CapPerCustomer = (ll-0.5);


end










