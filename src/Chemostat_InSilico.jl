module Chemostat_InSilico

    import ProjAssistant
    const PjAss = ProjAssistant
    PjAss.@gen_top_proj

    include("Dynamic/Dynamic.jl")

    function __init__()
        PjAss.@create_proj_dirs
    end

end
