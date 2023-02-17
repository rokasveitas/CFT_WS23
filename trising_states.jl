using Combinatorics, Serialization, Plots


let


states = []

# N is system size,
# ns is Vector of excitation momenta
function H(N, ns; pbc=false)
	# E0 + sum n_i ϵ_i
	
	E = 0.
	p = 0.
	ks = (2*π/N ).* [n for n in ns]

	if pbc
		E += -2 * (cot(π / (2*N)) - csc(π / (2*N)))
	else
		ks .-= π / N
	end


	for k in ks
		E += 4 * abs(sin(k/2))
		p += k
	end

	return E,p * N / (2*π)
end

N = 2000
pmax = 6
nmax = 3


for nparts = 0:2:nmax # number of particles
	for ns in Iterators.product([-pmax:pmax for part in 1:nparts]...)
		if unique(ns) == [n for n in ns]
			push!(states, (ns, H(N, ns; pbc=true), true))
			push!(states, (ns, H(N, ns; pbc=false), false))
		end
	end
end


sort!(states; by= x -> x[2][1])

println(length(states))


abcstates = filter(x -> !x[3], states);
pbcstates = filter(x -> x[3], states);

pbcground = pbcstates[1][2][1]

aens = [s[2][1] for s in abcstates]
pens = [s[2][1] for s in pbcstates]
amos = [s[2][2] for s in abcstates]
pmos = [s[2][2] for s in pbcstates]

aens /= pbcground * 8
pens /= pbcground * 8

scatter(amos, aens; color=:blue, label="ABC")
scatter!(pmos, pens; color=:red, label="PBC")

plot!(ylims=(0,5), xlims=(-.05,.05))

Plots.savefig("trising.pdf")

serialize("trising_states", states)

end