# Noah-MP version 4.3 release

This is an official Noah-MP source code for version 4.3, which is consistent with that released in WRF v4.3. 

## This new version includes the following major updates:

**(1) Snow-related update**:

Add wet-bulb temperature snow-rain partitioning scheme (OPT_SNF=5) based on Wang et al. 2019 (NWM)

Add snow retention process at the snowpack bottom to improve streamflow modeling (NWM)

Modify wind-canopy absorption coefficient (CWPVT) parameter values in MPTABLE to be vegetation dependent based on Goudriaan1977

Bring hard-coded snow emissivity and parameter (2.5*z0) in snow cover formulation to tunable MPTABLE parameters

Update MFSNO in snow cover formulation with optimized vegetation-dependent values

Limit the bulk leaf boundary layer resistance (RB) to a more realistic range (5~50)

**(2) New irrigation scheme**:

multiple irrigation methods: sprinkler, micro, and surface flooding

**(3) Crop scheme update**: 

separate the original generic crop physiology parameters in the modis vegetation section into C3/C4 specific parameters in the crop section

**(4) New urban physics working with Noah-MP**:

Local climate zone (LCZ), solar panel, green roof, new building drag parameterization
