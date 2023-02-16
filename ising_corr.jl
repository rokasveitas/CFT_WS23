using LinearAlgebra, Plots, Serialization, Combinatorics, PyCall

let

function uv(k)
	u = cos(k) - 1 - sqrt(2 - 2 * cos(k))
	v = sin(k)

	return k==0 ? [-1/sqrt(2), 1/sqrt(2)] : [u,v] / norm([u,v])
end


function prop(ks)
	k, q = ks

	return k == -q ? uv(q)[1]' * uv(-k)[2]' : 0.
end



pf = pyimport("pfapack.pfaffian")

function corr(ks)
	n = length(ks)

	props = zeros(n,n)

	for i=1:n
		for j = (i+1):n
			props[i,j] = prop([ks[i], ks[j]])
		end
	end

	@show prop([pi/2, -pi/2])
	@show props

	@show eigvals(props)

	props -= props'
	@show props


	@show eigvals(props)

	pf.pfaffian(props)
end


println(corr([0, pi, pi/2, -pi/2]))



end