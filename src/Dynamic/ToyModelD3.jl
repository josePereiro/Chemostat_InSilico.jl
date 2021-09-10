## -------------------------------------------------------------------
function ToyModel3D()
    net = Dict()

    # glyc
    # https://en.wikipedia.org/wiki/Glycolysis
    # Glucose + 2 NAD+ + 2 ADP + 2 Pi → 2 pyruvate + 2 NADH + 2 H+ + 2 ATP

    # ppp
    # https://en.wikipedia.org/wiki/Pentose_phosphate_pathway
    # The overall reaction for this process is: Glucose 6-phosphate + 2 NADP+ + H2O → ribulose 5-phosphate + 2 NADPH + 2 H+ + CO.

    # resp
    # https://en.wikipedia.org/wiki/Oxidative_phosphorylation
    # 1/2 O2 + NADH + H+ → H2O + NAD+
    # O2 + 2*NADH + 2*H+ → 2*H2O + 2*NAD+
    # When one NADH is oxidized through the electron transfer chain, three ATPs are produced
    # When one NAPDH is oxidized through the electron transfer chain, two ATPs are produced
    # we use 2.5 as average
    # O2 + 2*NADH + 2*H+ → 2*H2O + 2*NAD+ 5*ATP

    # tac
    # https://en.wikipedia.org/wiki/Citric_acid_cycle
    # Products of the first turn of the cycle are one GTP (or ATP), three NADH, one FADH2 and two CO2.
    # Because two acetyl-CoA molecules are produced from each glucose molecule, two cycles are required per glucose molecule.
    # Therefore, at the end of two cycles, the products are: two GTP, six NADH, two FADH2, and four CO2.

    # ferm
    # http://ecolistudentportal.org/article_fermentation#_
    # Many fermenting bacteria including E. coli can convert acetyl-CoA to 
    # acetate to generate ATP. The two key enzymes are phosphate acetyltransferase and 
    # acetate kinase. This mini-pathway generates 1 mole of acetate and 
    # 1 mole of ATP (by SLP) per acetyl-CoA consumed.
    # acetyl-CoA + Pi → acetyl-P + CoA
    # acetyl-P + ADP → acetate + ATP
    # AcCoa -> Ac + Atp

    # z biomass
    # Data From:
    # Varma, (1993): 2465–73. https://doi.org/10.1128/AEM.59.8.2465-2473.1993.
    # FIG 2
    # GBY = 0.068 # mmol GLC/gDW
    GBY = 1 / 0.068 # gDW/ mmol

    # feistGenomescaleMetabolicReconstruction2007 DOI: 10.1038/msb4100155
    # Using chemostat data for E. coli growing aerobically on
    # glucose (see Supplementary information), we estimated the
    # GAM and NGAM costs (Figure 3D). We found that an NGAM
    # value of 8.39 mmol ATP gDW1 h1 and a GAM value of
    # 59.81 mmol ATP gDW1 best fit the experimental data.
    GAM = 59.81 # EColi GAM 59.81 mmol ATP gDW1
    NGAM = 8.4 # EColi NGAM 8.4 mmolATP/gDW/h

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
    #= Atp   =#  +2.0    0.0   +5.0   +1.0   +1.0   -GAM  -NGAM   0.0    0.0    0.0 ;  #= Atp   =#
    #= Ac    =#   0.0    0.0    0.0    0.0   +1.0    0.0   0.0   +1.0    0.0    0.0 ;  #= Ac    =#
    #= NADH  =#  +4.0   +2.0   -2.0   +4.0    0.0    0.0   0.0    0.0    0.0    0.0 ;  #= NADH  =#
    #= AcCoa =#  +1.0   +1.0    0.0   -1.0   -1.0    0.0   0.0    0.0    0.0    0.0 ;  #= AcCoa =#
    #= Oxy   =#   0.0    0.0   -1.0    0.0    0.0    0.0   0.0    0.0    0.0   +1.0 ;  #= Oxy   =#
    ]
    
    net[:mets] = ["Glc", "Atp", "Ac", "NADH", "AcCoa", "Oxy"]
    net[:b] =    [0.0  , 0.0  , 0.0 , 0.0   , 0.0    , 0.0  ] # const exchanges
    
    
    net[:rxns] = ["glyc"  , "ppp"  , "resp" , "tac" , "ferm" , "z"   , "atpm", "ua"  ,  "ug" , "uo" ]
    net[:lb]   = [ 0.0    , 0.0    , 0.0    , 0.0   , 0.0    ,  0.0  ,  1.0  , -AB   ,  0.0  , 0.0  ]
    net[:ub]   = [ AB     , AB     , AB     , AB    , AB     ,  Vz   ,  1.0  ,  0.0  ,  Vg   ,  AB  ]
    net[:c]    = [ 0.0    , 0.0    , 0.0    , 0.0   , 0.0    ,  -1.0 ,  0.0  ,  0.0  ,  0.0  , 0.0  ]
    
    return reducebox!(MetNet(;net...))
end