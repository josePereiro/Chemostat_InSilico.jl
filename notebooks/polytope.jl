### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 8b576bd2-25d3-11eb-2a6d-3fc198380109
begin
	import Chemostat_Dynamics: SimpleDynamic
	SD = SimpleDynamic
	using Plots
	using PlutoUI
end

# ╔═╡ e763fdb4-25d3-11eb-1708-0dc36f56d335
begin 
	
	md"""
	Choose a set of parameters
	
	Lr $(@bind Lr Slider(-10.0:0.01:0.0; default=0, show_value=true))
	Ur $(@bind Ur Slider(0.0:0.01:10.0; default=0.45, show_value=true))
	
	Lg $(@bind Lg Slider(-10.0:0.01:0.0; default=0, show_value=true))
	Ug $(@bind Ug Slider(0.0:0.01:10.0; default=0.5, show_value=true))
	
	Ll $(@bind Ll Slider(-10.0:0.01:10.0; default=-10.0, show_value=true))
	Ul $(@bind Ul Slider(-10.0:0.01:10.0; default=0.1, show_value=true))
	
	Nbioms $(@bind Nb Slider(1.0:1000.0; default=348, show_value=true))
	Nf $(@bind Nf Slider(1.0:10.0; default=2, show_value=true))
	Nr $(@bind Nr Slider(1.0:50.0; default=18, show_value=true))
	
	cg $(@bind cg Slider(0.01:0.01:50.0; default=15, show_value=true))
	xi $(@bind xi Slider(0.001:0.001:500.0; default=50.0, show_value=true))
	"""
end

# ╔═╡ 38248926-25d4-11eb-3440-d5b8d471db91
SD.plot_polytop(SD.Polytope(;Lr ,Ur ,Lg ,Ug ,Ll ,Ul ,Nb ,Nf ,Nr ,cg ,xi))

# ╔═╡ dbc538ce-2608-11eb-16a1-79c5dd8b280d
x = 1

# ╔═╡ Cell order:
# ╟─e763fdb4-25d3-11eb-1708-0dc36f56d335
# ╟─38248926-25d4-11eb-3440-d5b8d471db91
# ╠═8b576bd2-25d3-11eb-2a6d-3fc198380109
# ╠═dbc538ce-2608-11eb-16a1-79c5dd8b280d
