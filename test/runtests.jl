using Chemostat_Dynamics
const SD = Chemostat_Dynamics.SimpleDynamic
using Test

@testset "Chemostat_Dynamics.jl" begin
    # Write your tests here.

    ## ------------------------------------------------------------------
    @testset "Polytope Bounds Tests" begin
        p = SD.SimParams() # Default params
        # The polytope is pointy in this regions
        @test SD.Δvg(SD.vatp_global_max(p), p) == 0.0
        @test SD.Δvg(SD.vatp_global_min(p), p) == 0.0
        @test SD.Δvatp(SD.vg_global_min(p), p) == 0.0

        # The polytope is NOT pointy in this regions_
        @test SD.Δvatp(SD.vg_global_max(p), p) > 0.0
    end

    ## ------------------------------------------------------------------
    @testset "Sampling Polytope Tests" begin
        n = Int(5e6)
        p = SD.SimParams() # Default params
        cells_pool = SD.generate_random_cells(p, n; tries = 100, verbose = true)
        @test all(SD.is_inpolytope.(cells_pool))
        @test all(map((cp) -> cp === p, getfield.(cells_pool, :p)))

        mvatp = SD.vatp_global_max(p)
        pcells = SD.pick_cells(n, cells_pool) do cell
            prob = SD.vatp(cell)/mvatp
            return rand() <= prob
        end
        @test all(SD.is_inpolytope.(pcells))
        @test all(map((cp) -> cp === p, getfield.(pcells, :p)))
        @test sum(SD.vatp.(cells_pool)) < sum(SD.vatp.(pcells)) # See picking function
    end

end


