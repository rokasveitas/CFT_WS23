using LinearAlgebra, Plots, Serialization, Combinatorics, PyCall


let


function uv(k)
	u = cos(k) - 1 - sqrt(2 - 2 * cos(k))
	v = sin(k)

	return k==0 ? [-1/sqrt(2), 1/sqrt(2)] : [u,v] / norm([u,v])
end

fourier_vec(k, N) = exp.(-im * k * (1:N)) / sqrt(N)

function γ_dec(k, N)
	U = zeros(Complex{Float64}, N)
	V = zeros(Complex{Float64}, N)
	
	if k == 0.
		fill!(U, 1/sqrt(N))
		return U, V
	end
	uvk = uv(k)

	U = fourier_vec(k, N) * uvk[2]
	V = conj(fourier_vec(-k,N)) * uvk[1]

	return U, V
end

function γdag_dec(k, N)
	U, V = γ_dec(k, N)

	return conj(V), conj(U)
end

pf = pyimport("pfapack.pfaffian")

function σσ(N)
	U = zeros(Complex{Float64}, 2*N, N)
	V = zeros(Complex{Float64}, 2*N, N)

	# <φ_I|
	for n=1:N÷2

		k = (n - 1/2) * 2*pi/N
		un, vn = γdag_dec(k, N)
		unegn, vnegn = γdag_dec(-k, N)

		U[2*n-1, :] = un
		U[2*n, :]   = unegn

		V[2*n-1, :] = vn
		V[2*n, :]   = vnegn
	end

	# c_1 + cdag_1
	U[N+1, 1] = 1.
	V[N+1, 1] = 1.

	# γ0
	γ0dec = γ_dec(0., N)
	U[N+2, :] = γ0dec[1]
	V[N+2, :] = γ0dec[2]

	# |φ_σ>
	for n=1:(N÷2 - 1)
		k = 2*pi*n / N
		unegn, vnegn = γ_dec(-k, N)
		un, vn = γ_dec(k, N)

		U[N+2+2*n-1, :] = unegn
		U[N+2+2*n  , :] = un

		V[N+2+2*n-1, :] = vnegn
		V[N+2+2*n  , :] = vn
	end

	A = UpperTriangular(U * V')
	A = A - transpose(A)
	# display(real.(A))
	# println()
	# @show norm(A)
	# @show norm(real.(A))
	# @show norm(imag.(A))
	# @show norm(A - A')
	# @show norm(A - transpose(A))
	# @show norm(A + A')
	# @show norm(A + transpose(A))
	display(eigvals(A))
	println()
	return pf.pfaffian(A) / sqrt(norm_pbc_bogvac(N)) / sqrt(norm_abc_bogvac(N))
end

function norm_pbc_bogvac(N)
	U = zeros(Complex{Float64}, N - 1, N)
	V = zeros(Complex{Float64}, N - 1, N)
	
	# γ0
	#γ0dec = γ_dec(0., N)
	#U[1, :] = γ0dec[1]
	#V[1, :] = γ0dec[2]

	for n=1:(N÷2 - 1)
		k = 2*pi*n / N
		unegn, vnegn = γ_dec(-k, N)
		un, vn = γ_dec(k, N)

		U[2*n  , :] = unegn
		U[2*n+1, :] = un

		V[2*n  , :] = vnegn
		V[2*n+1, :] = vn
	end

	fullU = vcat(conj(V)[end:-1:2, :], U[2:end,:])
	fullV = vcat(conj(U)[end:-1:2, :], V[2:end,:])


	A = UpperTriangular(fullU * fullV')
	A = A - transpose(A)
	println(pf.pfaffian(A))
	return pf.pfaffian(A)
end


function norm_abc_bogvac(N)
	U = zeros(Complex{Float64}, N , N)
	V = zeros(Complex{Float64}, N , N)
	
	# γ0
	#γ0dec = γ_dec(0., N)
	#U[1, :] = γ0dec[1]
	#V[1, :] = γ0dec[2]

	for n=1:(N÷2)
		k = 2*pi*(n - 1/2) / N
		unegn, vnegn = γ_dec(-k, N)
		un, vn = γ_dec(k, N)

		U[2*n-1, :] = unegn
		U[2*n  , :] = un

		V[2*n-1, :] = vnegn
		V[2*n  , :] = vn
	end

	fullU = vcat(conj(V)[end:-1:1, :], U[1:end,:])
	fullV = vcat(conj(U)[end:-1:1, :], V[1:end,:])


	A = UpperTriangular(fullU * fullV')
	A = A - transpose(A)
	println(pf.pfaffian(A))
	return pf.pfaffian(A)
end



@time println(σσ(4))


Ns = 4:2:20
σσs = σσ.(Ns)
println(σσs)

plot(abs.(σσs), Ns; xaxis=:log, yaxis=:log)



end