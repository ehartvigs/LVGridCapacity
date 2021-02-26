%% ------------------------------------------------------------------------
% This is the network generator, the inputs are population density, load
% profile for the customers, solar production profiles, and external model
% parameters found in the MainModelUK file.
%% ------------------------------------------------------------------------


function  [frac_AP,CustomersPerFeeder,MaxLoad,CustomersPerTransformer,TransformerType...
    ,LLVMax,solcap,NumberOfTransformers,CapPerCustomer,EnergyPerKM2,Limiter]...
    = NetworkModelUK(PopDensity, HouseLoadProfile, AptLoadProfile, SolProfile, Alpha, VoltageLimit)



%%  Set model parameters  %%
% Design voltage limits
Design_voltage = [1.1 0.94];   % Design voltage UK http://www.legislation.gov.uk/uksi/2002/2665/regulation/27/made

% Design loop impedance from IEC
Design_Z = 0.64*1000;

% Nominal voltage
Vn = 400;

% Demographic parameters
PeoplePerHouse = 2.35;
PeoplePerApt = 2.35;

% Fuse ratings (three phase equivalent)
AptFuseRating = 10;
HouseFuseRating = 20;

% Transformer parameters
TransformerCap = [1250 800 500 315 200 100];
TransformerCost = [195272 134751 101565 70501 53509 38446];
Transformer_Impedance = [6.5 10 13 20 32 65]/1000;
Z_TransformerR_list = [0.000848114120254821 0.00148841681507050 0.00251028652976876 0.00492770793724613 0.00776114000116266 0.0202385770250776 0.0404771540501553];
Z_TransformerX_list = [0.00763302708229339 0.0119073345205640 0.0125514326488438 0.0197108317489845 0.0310445600046506 0.0607157310752329 0.121431462150466];

% Cable parameters
CableCapacity = [52 67 80 94 138 178 230 345];
fuses_tripping_size = [6 10 16 20 25 32 40 50 63 80 100 125 160 200 250 315;...
    18 32 65 85 110 150 190 250 320 425 580 715 950 1250 1650 2200];
targetFuseRatings = [6 10 16 20 25 32 40 50 63 80 100 125 160 200 250 315];
Z_lineR_list = [1.83 1.15 0.76 0.641 0.32 0.206 0.125 0.062];       % Last value extrapolated
Z_lineX_list = [0.0861 0.0817 0.0783 0.0779 0.0762 0.0745 0.0752 0.0752];
Z_Line = sqrt(Z_lineR_list.^2+Z_lineX_list.^2);


% Tripping criteria parameters
c = 0.95;       
Ufn = 230;


% Set feeder length multiplier (gamma)
if PopDensity<=200
    LLAF_LV = 1; 
elseif PopDensity>200 && PopDensity<=1000
    LLAF_LV = 1.1;
else
    LLAF_LV = 1.2;
end


% Demand parameters
margin = Alpha;     % Transformer margin
AptBuilding_share = 1.25;     % Increased share of electricity consumption for
% apartment due to other appliances (e.g. elevator, lights and heating)
% https://www.sciencedirect.com/science/article/pii/S0378778811005019#fig0010
PowerFactor = 0.9;
PowerFactor = sind(acosd(PowerFactor));


% Cable costs per km (SEK)
if PopDensity>1000
    CostPerKmLineLV = 827000;
    CostPerKmLineMV = 1140746;
elseif PopDensity>200
    CostPerKmLineLV = 540000;
    CostPerKmLineMV = 887790;
else
    CostPerKmLineLV = 177000;
    CostPerKmLineMV = 512337;
end


% Calculate the share of single family and multifamily households.
if PopDensity <= 14690
    frac_AP_old = min([0.1173*log(0.34*PopDensity+47.8547) 1]);
    frac_AP = max([frac_AP_old 0]);
    frac_HH = 1-frac_AP_old;
else
    frac_AP = 1;
    frac_HH = 0;
end


fuse = AptFuseRating*frac_AP+HouseFuseRating*frac_HH;
PeoplePerHousehold = PeoplePerHouse*frac_HH + PeoplePerApt*frac_AP;


% Special cases of low population density
if PopDensity<4
    NumberOfCustomers = PopDensity;
    frac_AP = 0;
    frac_HH = 1;
else

    NumberOfCustomers = PopDensity/PeoplePerHousehold;
end



% Construct average load profile
AVG_LoadProfile = (frac_HH.*HouseLoadProfile+frac_AP.*AptBuilding_share*AptLoadProfile);  

maxP=fuse*Vn*sqrt(3)/1000;

% Set initial values for cost-minimization 
CustomersPerTransformer = NumberOfCustomers;
ADMDTr = 1.5*frac_AP+2.1*frac_HH;
PeakLoad = margin*0.7*CustomersPerTransformer*ADMDTr*(1+12/(ADMDTr*CustomersPerTransformer));
n0 = ceil(PeakLoad/TransformerCap(1));
MaxLoad = PeakLoad/margin;
ToTCost = zeros((200),length(TransformerCap));


% Finding lest cost option of transformer capacity and number of
% transformers to supply the load.
for gg = n0:200
    
    CustomersPerTransformer = round(NumberOfCustomers/gg);
    ADMDTr = 1.5*frac_AP+2.1*frac_HH; 
    PeakLoadTR = margin*0.7*CustomersPerTransformer*ADMDTr*(1+12/(ADMDTr*CustomersPerTransformer));
    indexTR = find((PeakLoadTR./(TransformerCap)<1)~=0,1, 'last');
    TransformerCostArea = TransformerCost(indexTR)*gg;
    
    A_TS = 1/gg; % in km^2
    NubmerOfConnectionsTransformer = NumberOfCustomers/gg;
        
    d = sqrt(A_TS)/(sqrt(NubmerOfConnectionsTransformer)+1);
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

    ToTCost((gg-n0+1),indexTR) = TransformerCostArea + LineCost; % Before 0.12
 
end

ToTCost(ToTCost==0)=NaN;

MinCost=min(min(ToTCost));
size(ToTCost);
[xt,TrSize]=find(ToTCost==MinCost);
xt=xt(1);
TrSize=TrSize(1);
NumberOfTransformers = n0+xt-1;


TransformerType = TransformerCap(TrSize);
A_TS = 1/NumberOfTransformers; % in km^2

CustomersPerTransformer= round(NumberOfCustomers/NumberOfTransformers);
ADMDTr =1.5*frac_AP+2.1*frac_HH;
PeakLoadTR = 0.7*CustomersPerTransformer*ADMDTr*(1+12/(ADMDTr*CustomersPerTransformer));
coincidenceTR = PeakLoadTR/(maxP*CustomersPerTransformer);


% Handle special cases of the uniform distribution of customers
if CustomersPerTransformer == 2
    LLVMax = 0.1;
    d = LLVMax;
    n = 2;
    nc = 2;
elseif CustomersPerTransformer == 1
    LLVMax =  0.100;
    d = LLVMax;
    n = 1;
    nc = 1;
else
    n = round(sqrt(CustomersPerTransformer));
    nc = round(CustomersPerTransformer/4);
    if n == 2
        d = sqrt(A_TS)/(4)+0.02;
        LLVMax = d;
        nc = 2;
    else
        LLVMax = LLAF_LV*sqrt(A_TS)*(n-1)/(n)+0.02; 
        d = LLVMax/n;

    end
end


% Calcualte transformer resistance and reactance
Z_transformer = Transformer_Impedance(TrSize);
Transformer_R = Z_TransformerR_list(TrSize);
Transformer_X = Z_TransformerX_list(TrSize);


% Set number of branches per feeder, max number of branches is 5
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
elseif n == 1
    m = zeros(1,1);
    mc = zeros(1,1);
    m(1) = n;
    mc(1) = nc;
    p = 1;
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
else
    m = zeros(1,2);
    mc = zeros(1,2);
    m(1)= round((n-1)/2);
    m(2) = n-m(1);
    mc(1)= round((nc-1)/2);
    mc(2) = nc-mc(1);
    p = 2;
end

% Calculate length of each feeder section, and aggregate for a total feeder
% length.
mc = flip(mc);
d_long = zeros(1,p);
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
R = zeros(1,p);
X = zeros(1,p);
ZmaxThick = zeros(9,p);
coincidenceLine = zeros(1,p);
P_demand = zeros(1,p);
Z_loop = [Z_transformer*1000];
L_max = zeros(9,p);
RX_multiplier = ones(1,p);
I_line_fuse = zeros(1,p);
mp = zeros(1,p);
z_supply = zeros(1,p+1);
targetLineImpedance = [4.18 2.63 1.91 1.47 0.746 0.495 0.324 0.324 0.324];

% ADMD, 1.5 and 2.1 are ADMD for apartments and houses respectively
ADMDLine = 1.5*frac_AP+2.1*frac_HH;
CustomersPerFeeder = sum(mc);


% Calculate initial feeder size
for rt=1:p
    rr = sum(mc(rt:end));
    mp(rt) = rr;
    
    PeakLoadLine = 0.7*rr*ADMDLine*(1+12/(ADMDLine*rr));
    coincidenceLine(rt) = PeakLoadLine/(rr*(0.7*ADMDLine*(1+12/(ADMDLine))));
    P_demand(rt) =  maxP*rr*coincidenceLine(rt);
    I_line = coincidenceLine(rt).*fuse*rr;
    
    % Check that cable current carrying capacity is higher than
    % cable fuse rating
    while I_line>max(targetFuseRatings)
        rr = round(rr/2);
        RX_multiplier(rt) = RX_multiplier(rt)*0.5;
        
        PeakLoadLine = 0.7*rr*ADMDLine*(1+12/(ADMDLine*rr));
        coincidenceLine(rt) = PeakLoadLine/(rr*(0.7*ADMDLine*(1+12/(ADMDLine))));
        P_demand(rt) =  maxP*rr*coincidenceLine(rt);
        I_line = coincidenceLine(rt).*fuse*rr;
    end
    
    
    % Check that earth fault impedance is small enough (and
    % that the length is not to long)
    I_line_fuse(rt) = targetFuseRatings(find((ceil(I_line)./targetFuseRatings<=1)~=0,1, 'first'));
    index = find((I_line_fuse(rt)./targetFuseRatings<=1)~=0,1, 'first') ;
    Iu = fuses_tripping_size(2,(index));
    ZmaxThick(rt) = 1000*c*Ufn/Iu;
    L_max(:,rt) = (ZmaxThick(rt)-sum(Z_loop))./targetLineImpedance';
    
    % Check that feeder length is within maximum allowed due to
    % tripping times
    while sum(L_max(:,rt)>d_long(rt)) == 0
        rr = round(rr/2);
        RX_multiplier(rt) = RX_multiplier(rt)*0.5;
        PeakLoadLine = 0.7*rr*ADMDLine*(1+12/(ADMDLine*rr));
        coincidenceLine(rt) = PeakLoadLine/(rr*(0.7*ADMDLine*(1+12/(ADMDLine))));
        P_demand(rt) =  maxP*rr*coincidenceLine(rt);
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
z_supply(1) = Z_transformer*1000;

% Check cable capacity and calcaulte supply impedance
for gg=1:p
    ixCable(gg) = find(L_max(:,gg),1);
    Cable(gg) = CableCapacity(ixCable(gg));
    if CableCapacity(ixCable(gg))<I_line_fuse(gg)
        ixCable(gg) = find((I_line_fuse(gg)./CableCapacity<=1)~=0,1, 'first');
        Cable(gg) = CableCapacity(ixCable(gg));
    end
    R(gg) = Z_lineR_list(ixCable(gg))*d_long(gg)*RX_multiplier(gg);
    X(gg) = Z_lineX_list(ixCable(gg))*d_long(gg)*RX_multiplier(gg);
    z_supply(gg+1) =  Z_Line(ixCable(gg))*d_long(gg)*RX_multiplier(gg);
end


id = find(ixCable==max(ixCable),1);

for hh=1:id
    Cable(hh) = CableCapacity(ixCable(id));
    R(hh) = Z_lineR_list(ixCable(id))*d_long(id);
    X(hh) = Z_lineX_list(ixCable(id))*d_long(id);
    z_supply(hh+1) = z_supply(hh) + Z_Line(ixCable(hh))*d_long(hh)*RX_multiplier(hh);
end

% Make sure updated resistance and reactance have correct quantity
R = R./1000;
X = X./1000;


% Calculate initial voltage drop
V = zeros(1,p+1);
V(1) = 400 -(coincidenceTR*CustomersPerTransformer*(Transformer_R.*(fuse*230*3)+Transformer_X.*(fuse*230*3).*PowerFactor)./Vn);
for kk=1:p
    V(kk+1) = V(kk)-(R(kk).*(P_demand(kk).*1000)+(X(kk)).*P_demand(kk).*1000*PowerFactor)./(Vn);
end

% Checking that voltage drop is within limits%, otherwise increse cable
% capacity
while  min(V./Vn)<Design_voltage(2)
    
    [val pos] = min(ixCable);
    if min(ixCable) == 8
        [val2 pos] = max(RX_multiplier);
        RX_multiplier(pos) = RX_multiplier(pos)*0.5;
        R(pos) = Z_lineR_list(val)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val;
        Cable(pos) = CableCapacity(val)*(1/RX_multiplier(pos));
        z_supply(pos+1) = z_supply(pos) + targetLineImpedance((val))*d_long(pos)*RX_multiplier(pos);
    else
        R(pos) = Z_lineR_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        X(pos) = Z_lineX_list(val+1)*d_long(pos)*RX_multiplier(pos)/1000;
        ixCable(pos) = val +1;
        Cable(pos) = CableCapacity(val+1);
        z_supply(pos+1) = z_supply(pos) + targetLineImpedance((val+1))*d_long(pos)*RX_multiplier(pos);
    end
    V = zeros(1,p+1);
    V(1) = 400 -min((coincidenceTR*CustomersPerTransformer*(Transformer_R.*(fuse*230*3)+Transformer_X.*(fuse*230*3).*PowerFactor)./Vn));
    for kk=1:p
        V(kk+1) = V(kk)-(R(kk).*(P_demand(kk).*1000)+(X(kk)).*P_demand(kk).*1000*PowerFactor)./(Vn);
    end
   
end


% Checking that supply impedance is below IEC values
while  sum(z_supply)>Design_Z
    
    [val pos] = min(ixCable);
    if min(ixCable) == 8
        [val2 pos] = max(RX_multiplier);
        RX_multiplier(pos) = RX_multiplier(pos)*0.5;
        z_supply(pos+1) =  Z_Line(val)*d_long(pos)*RX_multiplier(pos);
        ixCable(pos) = val;
        Cable(pos) = CableCapacity(val)*(1/RX_multiplier(pos));
    else
        z_supply(pos+1) =  Z_Line(val+1)*d_long(pos)*RX_multiplier(pos);
        ixCable(pos) = val+1;
        Cable(pos) = CableCapacity(val+1);
    end
    
    
end



% Add solar PV systems to each connection point on the feeder in increments
% of 0.5 kW
for ll=0.5:0.5:50
    
    % Calcualte maximum voltage along feeder
    V0 = 400 - min((CustomersPerTransformer*(Transformer_R.*(coincidenceTR*fuse*230*3-(ll).*SolProfile.*1000)+Transformer_X.*(coincidenceTR*fuse*230*3)*PowerFactor)./Vn));
    if length(mp) == 1
        voltage = V0- min((mp.*(R.*(coincidenceLine.*AVG_LoadProfile-(ll).*SolProfile.*1000)+X.*(coincidenceLine.*AVG_LoadProfile*PowerFactor))./Vn)');
    else
        voltage = V0- min(sum((mp.*(R.*(coincidenceLine.*AVG_LoadProfile-(ll).*SolProfile.*1000)+X.*(coincidenceLine.*AVG_LoadProfile*PowerFactor))./Vn)'));
    end
    
    % Calculate maximum demand along feeder
    deltaPowerPerCustomerDouble = mc.*max(abs(coincidenceLine(1).*AVG_LoadProfile-ll.*SolProfile.*1000))/400/sqrt(3);
    deltaPowerPerCustomer = max(abs(coincidenceTR*AVG_LoadProfile-ll.*SolProfile.*1000));
    deltaPower = deltaPowerPerCustomer.*CustomersPerTransformer;
    
    
    if (max(max(voltage))/400)>VoltageLimit(1)
        solcap = (ll-0.5)*NumberOfCustomers;
        Limiter = 10;
        break

    elseif deltaPowerPerCustomerDouble>Cable 
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










