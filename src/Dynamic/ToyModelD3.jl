## -------------------------------------------------------------------
function ToyModel3D()
    net = Dict()

    # feistGenomescaleMetabolicReconstruction2007 DOI: 10.1038/msb4100155
    # Using chemostat data for E. coli growing aerobically on
    # glucose (see Supplementary information), we estimated the
    # GAM and NGAM costs (Figure 3D). We found that an NGAM
    # value of 8.39 mmol ATP gDW1 h1 and a GAM value of
    # 59.81 mmol ATP gDW1 best fit the experimental data.
    GAM = 59.81 # EColi GAM 59.81 mmol ATP gDW1
    NGAM = 8.4 # EColi NGAM 8.4 mmolATP/gDW/h

    # https://en.wikipedia.org/wiki/Glycolysis
    # The overall process of glycolysis is:
    # Glucose + 2 NAD+ + 2 ADP + 2 Pi → 2 pyruvate + 2 NADH + 2 H+ + 2 ATP
    # Glucose + 2 NAD+ + 2 ADP + 2 Pi → 2 ac-Coa + 2 NADH + 2 H+ + 2 ATP
    GEY = 2.0 # (glycolisis energy yield) Per glucose molecule

    # https://en.wikipedia.org/wiki/Citric_acid_cycle
    # The theoretical maximum yield of ATP through oxidation of one molecule of glucose in 
    # glycolysis, citric acid cycle, and oxidative phosphorylation is 38 (assuming 3 molar 
    # equivalents of ATP per equivalent NADH and 2 ATP per UQH2). 
    REY = (38.0 - GEY)/2 # (Respiratory energy yield) Per acetate molecule
    # The P/O ration is tuned to make the oxygen uptake compatible with
    # Varma, (1993): 2465–73. https://doi.org/10.1128/AEM.59.8.2465-2473.1993.
    # FIG 2
    RNO = 6.0

    # http://ecolistudentportal.org/article_fermentation#_
    # Many fermenting bacteria including E. coli can convert acetyl-CoA to 
    # acetate to generate ATP. The two key enzymes are phosphate acetyltransferase and 
    # acetate kinase. This mini-pathway generates 1 mole of acetate and 
    # 1 mole of ATP (by SLP) per acetyl-CoA consumed.
    # ferm: acetyl-CoA + Pi → acetyl-P + CoA
    # ferm: acetyl-P + ADP → acetate + ATP
    # ferm: AcCoa -> Ac + Atp

    # z biomass
    # Data From:
    # Varma, (1993): 2465–73. https://doi.org/10.1128/AEM.59.8.2465-2473.1993.
    # FIG 2
    # GBY = 0.068 # mmol GLC/gDW
    GBY = 1 / 0.068 # gDW/ mmol

    # uptakes limits
    AB = 1000.0 # abs bound
    # This model is bounded by the maximum rates found for EColi.
    # Data From:
    # Varma, (1993): 2465–73. https://doi.org/10.1128/AEM.59.8.2465-2473.1993.
    # Extract max exchages from FIG 3 to form the maximum polytope
    Vg = 20      # (glucose uptake limit) mmol/gDW h
    Vz = 0.5     # D < 0.5 below the acetate switch [23] Nature, 528(7580):99–104, 2015.

    net[:S] = [      
    #             glyc   ppp    resp   tac    ferm   z     atpm   ua     ug     uo    
    #= Glc   =#  -1.0   -1.0    0.0    0.0    0.0   -GBY   0.0    0.0   +1.0    0.0 ;  #= Glc   =#
    #= Atp   =#  +1.0    0.0   +1.0   +1.0    0.0   -GAM  -1.0    0.0    0.0    0.0 ;  #= Atp   =#
    #= Ac    =#   0.0    0.0    0.0    0.0   +1.0    0.0   0.0   +1.0    0.0    0.0 ;  #= Ac    =#
    #= NADH  =#  +1.0   +1.0   -1.0   +4.0    0.0    0.0   0.0    0.0    0.0    0.0 ;  #= NADH  =#
    #= AcCoa =#  +1.0   +1.0    0.0   -1.0   -1.0    0.0   0.0    0.0    0.0    0.0 ;  #= AcCoa =#
    #= Oxy   =#   0.0    0.0   -1.0    0.0    0.0    0.0   0.0    0.0    0.0   +1.0 ;  #= Oxy   =#
    ]
    
    net[:mets] = ["Glc", "Atp", "Ac", "NADH", "AcCoa", "Oxy"]
    net[:b] =    [0.0  , 0.0  , 0.0 , 0.0   , 0.0    , 0.0  ] # const exchanges
    
    
    net[:rxns] = ["glyc"  , "ppp"  , "resp" , "tac" , "ferm" , "z"   , "atpm", "ua"  ,  "ug" , "uo" ]
    net[:lb]   = [ 0.0    , 0.0    , 0.0    , 0.0   , 0.0    ,  0.0  , NGAM  , -AB   ,  0.0  , 0.0  ]
    net[:ub]   = [ AB     , AB     , AB     , AB    , AB     ,  Vz   , NGAM  ,  0.0  ,  Vg   ,  AB  ]
    net[:c]    = [ 0.0    , 0.0    , 0.0    , 0.0   , 0.0    ,  -1.0 , 0.0   ,  0.0  ,  0.0  , 0.0  ]
    
    return reducebox!(MetNet(;net...))
end