using LinearAlgebra, Plots, Serialization, Combinatorics, PyCall


let

# γ_k = uk c_k + vk cdag_-k
function uv(k; h=0)
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

function Base.:*(g::FermiOp, a::Complex)
	return FermiOp(g.U * a, g.V * a, g.N)
end

function Base.:/(g::FermiOp, a::Complex)
	return FermiOp(g.U / a, g.V / a, g.N)
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

	@show norm2
	return sqrt(real(norm2))
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


function equaltime_g(N, j)

	num = four_σ(N, j)

	den = σσ(N)^2

	return num / den
end

@time println(γsnorm(bogvacstring(4, false)))
@time println(γsnorm(bogvacstring(4, true)))

@time println(σσ(4))


Ns = 4:2:20
σσs = σσ.(Ns)
println(σσs)

plot(abs.(σσs), Ns; xaxis=:log, yaxis=:log)



end