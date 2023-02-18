include("IsingRadQuant.jl")

using .IsingRadQuant
using Statistics, Plots

let 

# Simple regressions to extract data
# y = m*x + b
# returns (m, b)
function linreg(x, y)
	m = sum((x .- mean(x)) .* (y .- mean(y))) / sum((x .- mean(x)) .^2)
	return m, mean(y) - m * mean(x)
end

# y = A * x^Δ
# returns (Δ, A)

# log y = log A + D log x
function powerreg(x, y)
	m, b = linreg(log.(x), log.(abs.(y)))
	return m, exp(b)
end



# First calculate <0|σ_1^x|φ_σ> ~ A / N^(Δ_σ) ∝ <σσ> to find that Δ_σ = 1/8.

Ns = 100:4:148

@time σσs = σσ.(Ns) # 6.490011 sec, 47.37% compilation

#scatter(Ns, abs.(σσs); yaxis=:log, xaxis=:log)
#savefig("ising_σσ.pdf")

Δσ, A =  powerreg(Ns, σσs)
@show Δσ #  -0.12499824248823803  ≈ -1/8
@show A  #   1.0105297694978246     (non-universal)


# Next calculate <φ_σ|σ_1^x|φ_ε> ~ A * <σσε> = A λ_σσε / N^(Δ_σ) to find λ_σσε = 1/2

@time λ = abs(σσε(100) / σσ(100)) # 0.393167 sec, 47.49% compilation

@show λ  #  0.49996915621782684 ≈ 1/2


# Now calculate <φ_σ|σ_1^x σ_j^x|φ_σ> ~ <σσσσ> ~ A^2 g(z,zbar) / N^(2Δ), so g(z,zbar) = <σσσσ> / <σσ>^2.  This only works on unit circle, z*zbar=1

N = 100

@time gs_num = [equaltime_g(N, j) for j=1:100] # 31.949648 sec, 21.05% compilation

θs = 0:0.001:2*pi
@time gs_ana = [g_polar(1., θ) for θ in θs] # 0.162349 sec, 97.27% compilation

@time scatter((1:100)*2*pi/100 .- 2*pi/100, abs.(gs_num); label="Numerical")
@time plot!(θs, abs.(gs_ana); color=:red, label="Exact")

@time savefig("equaltime_gs_test.pdf") # This plot shows both the numerical and exact results superimposed

scatter((1:100)*2*pi/100 .- 2*pi/100,abs.(gs_num) .- abs.([g_polar(1., θ) for θ in (1:100)*2*pi/100 .- 2*pi/100])) 
savefig("equaltime_gs_err.pdf") # This plot shows the difference in magnitudes between the numerical and exact results






end