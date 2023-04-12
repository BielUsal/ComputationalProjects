S1(G) = 10*(exp(im*pi*(G[1]+G[2]))+exp(im*pi*(G[2]+G[3]))+exp(im*pi*(G[1]+G[3]))) +10 *(exp(im*pi*G[1])+exp(im*pi*G[2])+exp(im*pi*G[3]))
S2(G)= 18*(exp(im*pi*(G[1]+G[2]))+exp(im*pi*(G[2]+G[3]))+exp(im*pi*(G[1]+G[3]))) +10*(exp(im*pi*G[1])+exp(im*pi*G[2])+exp(im*pi*G[3]))

Gengus = [[1,1,1],[2,0,0], [2,2,0],[3,1,1],[4,0,0],[3,3,1],[4,2,0]]
I1 = abs.(S1.(Gengus))
I2 = abs.(S2.(Gengus))

I1[2] *= 3
I1[3] *=2
I1[4] *=2
I1[5] *=3
I1[6] *= 3
I1[7] *=3
I2[2] *= 3
I2[3] *=2
I2[4] *=2
I2[5] *=3
I2[6] *= 3
I2[7] *=3
I1 .*= 3/maximum(I1)
I2 .*= 3/maximum(I2)

println("NaF:",I1)
println("NaCl", I2)
scatter(I1,normalize=:pdf)
scatter!(I2,normalize=:pdf)