using Plots
using Roots
β(n) = 11-2*n/3  
αₛ(alpha0,mu0,mu,n0,n) = 2*pi*alpha0/(2*pi+alpha0*(β(n)*log(mu)-β(n0)*log(mu0)))
α₀ = 0.1181
mz = 91.2 #GeV
mb = 4.2 #GeV
mc = 1.3 #GeV
nz = 5
nb = 4
nc = 3
α₁= αₛ(α₀,mz,mb,nz,nb)
println("α₁ = ",α₁)
α₂= αₛ(α₁,mb,mc,nb,nc)
println(α₂)
α₃(μ)= αₛ(α₂,mc,μ,nc,3)
plot(α₃,0.12666,0.5,label="α₃(μ)",xlabel="μ(GeV)",ylabel="α₃(μ)")
Λ = find_zero(μ->1/α₃(μ),0.1)
vline!([Λ],label="Λ = $Λ GeV")
savefig("renorm.png")