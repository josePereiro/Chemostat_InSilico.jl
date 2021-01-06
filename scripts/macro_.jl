using ProgressMeter

## ----------------------------------------------------------------------------
let
    @time begin
        N = 10
        prog = Progress(N; dt = 1.0);
        for i = 1:N

            # Computation
            sleep(3/N)
            next!(prog; showvalues = [
                    ("heavy val", (sleep(3/N); i)) 
                ]
            )
        end
    end
    finish!(prog)
end

## ----------------------------------------------------------------------------
let
    prog = Progress(1);
    @edit  next!(prog)
end

## ----------------------------------------------------------------------------
#=
Hi,

Sometimes I need to show data that is computed on the fly, because `next!` and `update!` are functions theirs arguments are evaluated every time, even in a lazy iteration, so the computation of `showvalues` is done

Example 
```julia
@time begin
        N = 10
        prog = Progress(N; dt = 1.0);
        for i = 1:N

            # Computation
            sleep(3/N)
            next!(prog; showvalues = [
                    ("heavy val", (sleep(3/N); i)) 
                ]
            )
        end
end
finish!(prog)
# Progress: 100%|███████████████████████| Time: 0:00:06
# heavy val:  10
# 6.133779 seconds (539 allocations: 38.859 KiB)
```
=#
## ----------------------------------------------------------------------------
macro mytime(ex)
    quote
        local elapsedtime = time_ns()
        local val = $(esc(ex))
        elapsedtime = time_ns() - elapsedtime
        println(elapsedtime)
        val
    end
end

## ----------------------------------------------------------------------------
macro test(ex)
    sleep(x) = sleep(x)
    quote
        $(esc(ex))
    end
end
## ----------------------------------------------------------------------------
@test sleep(1)
## ----------------------------------------------------------------------------
@macroexpand @test sleep(1)
## ----------------------------------------------------------------------------

