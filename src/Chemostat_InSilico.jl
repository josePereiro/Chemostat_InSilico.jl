module Chemostat_InSilico

    import UtilsJL
    const Ass = UtilsJL.ProjAssistant
    Ass.@gen_top_proj

    include("Dynamic/Dynamic.jl")

    function __init__()
        Ass.@create_proj_dirs
    end

end
