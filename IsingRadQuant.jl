

module IsingRadQuant

using LinearAlgebra, Plots, Serialization, Combinatorics, PyCall

export uv, FermiOp, γ, γdag, γstring, bogvacstring, γsnorm, timeevstr
export σσ, ε_state, σσε, two_σ, four_σ, equaltime_g, g_full, g_exact, g_polar

# bogolyubov transformation
# γ_k = uk c_k + vk cdag_-k
function uv(k)
	v = cos(k) - 1  - sqrt(2 - 2 * cos(k))
	u = sin(k)

	return k==0 ? [1/sqrt(2), -1/sqrt(2)] : [u,v] / norm([u,v])
end

struct FermiOp
	U::Vector
	V::Vector
	N::Integer
	
end

function γ(k, N)
	U = zeros(Complex{Float64}, N)
	V = zeros(Complex{Float64}, N)
	
	if k == 0.
		fill!(V, 1/sqrt(N))
		return FermiOp(U, V, N)
	end
	uvk = uv(k)

	U = fourier_vec(k, N) * uvk[2]
	V = conj(fourier_vec(-k,N)) * uvk[1]

	return FermiOp(U, V, N)
end

function γdag(k, N)
	γdagop = γ(k, N)
	return FermiOp(conj(γdagop.V), conj(γdagop.U), N)
end

fourier_vec(k, N) = exp.(-im * k * (1:N)) / sqrt(N)


pf = pyimport("pfapack.pfaffian")

function γstring(γs)
	U = vcat([transpose(g.U) for g in γs]...)
	V = vcat([transpose(g.V) for g in γs]...)

	A = UpperTriangular(U * transpose(V))
	A = A - transpose(A)

	return pf.pfaffian(A)
end

function γdag(γs::Vector)
	newγs = reverse(γs)
	for (i,g) in enumerate(newγs)
		newγs[i] = FermiOp(conj(g.V), conj(g.U), g.N)
	end

	return newγs
end

function Base.:*(g::FermiOp, a)
	return FermiOp(g.U * a, g.V * a, g.N)
end

function Base.:*(a, g::FermiOp)
	return FermiOp(g.U * a, g.V * a, g.N)
end

function Base.:/(g::FermiOp, a)
	return FermiOp(g.U / a, g.V / a, g.N)
end

function Base.:+(g::FermiOp, h::FermiOp)
	@assert g.N == h.N
	return FermiOp(g.U + h.U, g.V + h.V, g.N)
end

function Base.:-(g::FermiOp, h::FermiOp)
	@assert g.N == h.N
	return FermiOp(g.U .- h.U, g.V .- h.V, g.N)
end

γdag(g::FermiOp) = FermiOp(conj(g.V), conj(g.U), g.N)

function bogvacstring(N, pbc=true)
	if pbc
		γs = FermiOp[ γ(0., N) ]

		for n in 1:(N÷2-1)
			push!(γs, γ(-2*pi*n/N, N) / uv(2*pi*n/N)[2])
			push!(γs, γ( 2*pi*n/N, N) / uv(2*pi*n/N)[1])
		end
		return γs
	else
		γs = FermiOp[]

		for n in 1:N÷2
			push!(γs, γ(-(n - 1/2) * 2*pi/N, N) / uv((n - 1/2) * 2*pi/N)[2])
			push!(γs, γ( (n - 1/2) * 2*pi/N, N) / uv((n - 1/2) * 2*pi/N)[1])
		end

		return γs
	end
end

function γsnorm(γs)
	norm2 = γstring(vcat(γdag(γs), γs))

	#@show norm2
	return sqrt(real(norm2))
end


# Returns the string over the bc basis that evolves by imag. time τ=t
function timeevstr(N, t; pbc=false)
	ks = [(n - (pbc ? 0 : 1/2))*2*pi/N for n in 1:N]
	
	Astring = [(exp(t * 4 * abs(sin(k/2))) - 2) * γdag(γ(k, N)) + γ(k, N) for k in ks]
	Bstring = [exp(-t/2 * 4 * abs(sin(k/2))) * (γdag(γ(k, N)) - γ(k, N))  for k in ks]

	return [Astring Bstring][:]
end

function σσ(N)
	γs = γdag(bogvacstring(N, false))
	# σ1 = c_1 + cdag_1
	σ1 = FermiOp([1., zeros(N-1)...], [-1. * im, zeros(N-1)...], N)
	push!(γs, σ1)
	append!(γs, bogvacstring(N, true))

	return γstring(γs) / γsnorm(bogvacstring(N, false)) / γsnorm(bogvacstring(N, true))
end


function ε_state(N)
	γs = bogvacstring(N, false)

	insert!(γs, 1, γdag(γ(-pi/N, N)))
	insert!(γs, 1, γdag(γ( pi/N, N)))

	return γs
end

function σσε(N)
	γs = γdag(bogvacstring(N, true))

	σ1 = FermiOp([im, zeros(N-1)...], [1., zeros(N-1)...], N)
	push!(γs, σ1)

	append!(γs, ε_state(N))

	return γstring(γs) / γsnorm(bogvacstring(N, true)) / γsnorm(ε_state(N))
end

function two_σ(N, j)
	γs = γdag(bogvacstring(N, false))

	σ1 = FermiOp([im, zeros(N-1)...], [1., zeros(N-1)...], N)
	σj = FermiOp([zeros(j-1)..., im, zeros(N-j)...], [zeros(j-1)..., 1., zeros(N-j)...], N)
	push!(γs, σ1, σj)

	Astring = [FermiOp([zeros(k-1)..., im, zeros(N-k)...], [zeros(k-1)..., 1., zeros(N-k)...], N) for k=1:j-1]
	Bstring = [FermiOp([zeros(k-1)..., im, zeros(N-k)...], [zeros(k-1)...,-1., zeros(N-k)...], N) for k=1:j-1]

	σzstring = [Astring Bstring][:]
	append!(γs, σzstring)

	append!(γs, bogvacstring(N, false))

	return γstring(γs) / γsnorm(bogvacstring(N, false))^2
end


function four_σ(N, j)
	γs = γdag(bogvacstring(N, true))

	σ1 = FermiOp([im, zeros(N-1)...], [1., zeros(N-1)...], N)
	σj = FermiOp([zeros(j-1)..., im, zeros(N-j)...], [zeros(j-1)..., 1., zeros(N-j)...], N)
	push!(γs, σ1, σj)

	Astring = [FermiOp([zeros(k-1)..., im, zeros(N-k)...], [zeros(k-1)..., 1., zeros(N-k)...], N) for k=1:j-1]
	Bstring = [FermiOp([zeros(k-1)..., im, zeros(N-k)...], [zeros(k-1)...,-1., zeros(N-k)...], N) for k=1:j-1]

	σzstring = [Astring Bstring][:]
	append!(γs, σzstring)

	append!(γs, bogvacstring(N, true))

	return γstring(γs) / γsnorm(bogvacstring(N, true))^2
end

function four_σ(N, j, t)
	γs = γdag(bogvacstring(N, true))

	append!(γs, γdag(timeevstr(N, -t, pbc=false)))
	σ1 = FermiOp([im, zeros(N-1)...], [1., zeros(N-1)...], N)
	append!(γs, timeevstr(N, t, pbc=false))

	σj = FermiOp([zeros(j-1)..., im, zeros(N-j)...], [zeros(j-1)..., 1., zeros(N-j)...], N)
	push!(γs, σ1, σj)

	Astring = [FermiOp([zeros(k-1)..., im, zeros(N-k)...], [zeros(k-1)..., 1., zeros(N-k)...], N) for k=1:j-1]
	Bstring = [FermiOp([zeros(k-1)..., im, zeros(N-k)...], [zeros(k-1)...,-1., zeros(N-k)...], N) for k=1:j-1]

	σzstring = [Astring Bstring][:]
	append!(γs, σzstring)

	append!(γs, bogvacstring(N, true))

	return γstring(γs) / γsnorm(bogvacstring(N, true))^2

end


function equaltime_g(N, j)

	num = four_σ(N, j)

	den = σσ(N)^2

	return num / den
end

function g_full(N, j, t)
	num = four_σ(N, j, t)

	den = σσ(N)^2

	return num / den
end

function g_exact(Z, Zb)
	Z, Zb = Z + 0.0im, Zb +0.0im
	out  = sqrt(1. - sqrt(1. - Z)) * sqrt(1. - sqrt(1. - Zb)) / 2 / (1. - Z)^(1/8) / (1. - Zb)^(1/8)
	out += sqrt(1. + sqrt(1. - Z)) * sqrt(1. + sqrt(1. - Zb)) / 2 / (1. - Z)^(1/8) / (1. - Zb)^(1/8)
	return out
end

g_polar(r, θ) = g_exact(r * exp(im * θ), r * exp(-im * θ))


end