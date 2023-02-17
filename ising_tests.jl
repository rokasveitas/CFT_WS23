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

@time σσs = σσ.(Ns)

#scatter(Ns, abs.(σσs); yaxis=:log, xaxis=:log)
#savefig("ising_σσ.pdf")

Δσ, A =  powerreg(Ns, σσs)
@show Δσ #  -0.12499824248823803  ≈ -1/8
@show A  #   1.0105297694978246     (non-universal)


# Next calculate <φ_σ|σ_1^x|φ_ε> ~ A * <σσε> = A λ_σσε / N^(Δ_σ) to find λ_σσε = 1/2

@time λ = abs(σσε(100) / σσ(100))

@show λ 








end