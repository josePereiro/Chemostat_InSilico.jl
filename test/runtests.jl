using Chemostat_Dynamics
using Chemostat_Dynamics.Utilities
using Chemostat_Dynamics.Polytopes
using Chemostat_Dynamics.MonteCarlo
using Chemostat_Dynamics.MaxEnt
using Test

@testset "Chemostat_Dynamics.jl" begin
    # Write your tests here.

    ## ------------------------------------------------------------------
    @testset "Polytope" begin
        p = Polytope() # Default params
        # The polytope is pointy in this regions
        @test Δvg(vatpU(p), p) == 0.0
        @test Δvg(vatpL(p), p) == 0.0
        @test Δvatp(vgL(p), p) == 0.0

        # The polytope is NOT pointy in this regions_
        @test Δvatp(vgU(p), p) > 0.0
    end

    ## ------------------------------------------------------------------
    @testset "Monte Carlos" begin
        n = Int(5e6)
        p = Polytope() # Default params
        cells_pool = generate_random_cells(p, n; tries = 100, verbose = true)
        @test all(is_inpolytope.(cells_pool))
        @test all(map((cp) -> cp === p, getfield.(cells_pool, :p)))

        mvatp = vatpU(p)
        pcells = pick_cells(n, cells_pool) do cell
            prob = vatp(cell)/mvatp
            return rand() <= prob
        end
        @test all(is_inpolytope.(pcells))
        @test all(map((cp) -> cp === p, getfield.(pcells, :p)))
        @test sum(vatp.(cells_pool)) < sum(vatp.(pcells)) # See picking function
    end

    ## ------------------------------------------------------------------
    @testset "Utilities" begin
        @testset "Get chunck " begin
            for it in 1:100
                r = 1:rand(50:100)
                n = rand(3:5)
                chuncks = get_chuncks(r, n; th = 0)
                @test all(vcat(chuncks...) .== collect(r))
            end
        end
    end

    ## ------------------------------------------------------------------
    @testset "MaxEnt" begin
        @testset "vatp_marginal_pdf " begin
            for it in 1:100
                for marginal_pdf in [vatp_marginal_pdf, vg_marginal_pdf]
                    xi = rand()*1e3 + 3e-2
                    pol = Polytope(;xi)
                    beta = rand()*1e3 + 1e-3
                    rang, probs = marginal_pdf(pol, beta; n = Int(1e5))
                    @test length(rang) == length(probs)
                    @test isapprox(sum(probs .* step(rang)), 1.0; atol = 1e-7)
                end
            end
        end
    end
end


