using LinearAlgebra
using Unitful
using UnitfulAstro
using OrdinaryDiffEq
using Plots
using FLoops
using Distributed
#include("./GalaxyCollisionFunctions.jl")
#addprocs(4)
theme(:dracula)
const G = 4.3009E-3 *1u"pc *(km/s)^2 /Msun"
function format_parameters(galaxy_args)
    #I didn't have this at first, but this makes all galaxies have an uniform format. Not too crazy
    return Dict(
        "mass" => galaxy_args[1]*1u"Msun",
        "radius"     => galaxy_args[2]*1u"kpc",
        "center_pos" => galaxy_args[3].*1u"kpc",
        "center_vel" => galaxy_args[4].*1u"km/s",
        "normal"     => galaxy_args[5],
        "Nᵣ"    => galaxy_args[6],
        "Nₛ"    => galaxy_args[7],
        "softening"  => galaxy_args[8]
    )
end
function init_disk!(galaxy,dT=1E-4u"yr")
    #=
    This function takes a 'galaxy' as an argument, which is an array of 8 arguments=>
        [mass,radius,center_pos,center_vel,normal,N_rings,N_stars,softening]
    And outputs the star positions, velocities and the velocity scale.
    =#
    dr = (1-galaxy["softening"]) * galaxy["radius"]/galaxy["Nᵣ"] # Disregarding the softening, this is just the width of the rings(raidus/number of slices)
    StarsPerRing = trunc(Int64,galaxy["Nₛ"]/galaxy["Nᵣ"])
#-------------------------------------------------------------Rotations--------------------------------------------------------------
    if norm(galaxy["normal"]) == 0 
        Rotation = I
    else 
        cosθ = normalize(galaxy["normal"])[3]
        sinθ = √(1-cosθ^2)
        u = [0,0,1] × normalize(galaxy["normal"])
        if norm(u) == 0 
            Rotation = I
        else 
         u = normalize(u)

        Rotation = [
            u[1]*u[1]*(1-cosθ)+cosθ u[1]*u[2]*(1-cosθ)-u[2]*sinθ u[1]*u[3]*(1-cosθ)+u[1]*sinθ;

            u[2]*u[1]*(1-cosθ)+u[3]*sinθ u[2]*u[2]*(1-cosθ)+cosθ u[2]*u[3]*(1-cosθ)-u[1]*sinθ;

            u[3]*u[1]*(1-cosθ)+u[2]*sinθ u[3]*u[1]*(1-cosθ)+u[1]*sinθ u[3]*u[3]*(1-cosθ)+cosθ
            ]
        end
    end
#-----------------------------------------------------------------------------------------------------------------------------------
    galaxy["star_pos"] = []
    galaxy["star_vel"] = []

    R = galaxy["softening"] * galaxy["radius"]
    for i ∈ 1:galaxy["Nᵣ"] #iterate over the rings
        #Randomly distribute the initial radii and angles 
        rₛ =(R * ones(StarsPerRing) + dr *rand(StarsPerRing))./1u"pc" #creates an n-dimensional vector of radii from R to R+dr 
        ϕₛ = 2π *rand(StarsPerRing)

        # Positions 
        vecᵣ = (Rotation * ([rₛ.*cos.(ϕₛ),rₛ.*sin.(ϕₛ),zeros(StarsPerRing)])).*1u"pc"
        x = ustrip.(u"m",galaxy["center_pos"][1].+vecᵣ[1])#here, we strip our units into our 'integrating' unitlessness because units don't really play nice in arrays of arrays
        y = ustrip.(u"m",galaxy["center_pos"][2].+vecᵣ[2])
        z = ustrip.(u"m",galaxy["center_pos"][3].+vecᵣ[3])

        # Velocities
        Tₛ = 2π * uconvert.(u"s",sqrt.((rₛ*1u"pc").^3/(G *Introoder["mass"])))

        Δϕ = 2π *uconvert(u"s",dT)./Tₛ 

        vecᵥ = (Rotation* [(rₛ/(uconvert(u"s",dT)*1u"1/s")).*(cos.(ϕₛ)-cos.(ϕₛ-Δϕ)),rₛ/(uconvert(u"s",dT)*1u"1/s").*(sin.(ϕₛ)-sin.(ϕₛ-Δϕ)),zeros(StarsPerRing)])*1u"pc/s"
        v₁ = ustrip.(u"m/s",galaxy["center_vel"][1].+vecᵥ[1])
        v₂ = ustrip.(u"m/s",galaxy["center_vel"][2].+vecᵥ[2])
        v₃ = ustrip.(u"m/s",galaxy["center_vel"][3].+vecᵥ[3])
        for j ∈ 1:StarsPerRing
            push!(galaxy["star_pos"],[x[j],y[j],z[j]])
            push!(galaxy["star_vel"],[v₁[j],v₂[j],v₃[j]])
        end
        R += dr
    end 
    galaxy["vel_scale"]=uconvert(u"km/s",√(G*galaxy["mass"]/(0.5*R)))
end


function evolve_disk(galaxy,dT=1e-4u"yr",N_steps=100000,frames=500)
    divs = trunc(Int64,N_steps/frames)
    #Integration stuff
    dT=ustrip(u"s",dT)
    rₘ=ustrip(u"m",galaxy["softening"]*galaxy["radius"])
    Nₛ=galaxy["Nₛ"]
    
    #Galaxy stuff
    M = ustrip(u"Msun",galaxy["mass"])
    R₀ = ustrip.(u"m",galaxy["center_pos"])
    V₀ = ustrip.(u"m/s",galaxy["center_vel"])

    #star stuff
    rₛ=galaxy["star_pos"]
    vₛ=galaxy["star_vel"]

    function Gravity(dq::Vector{Float64},q::Vector{Float64},m,t::Float64)::Vector{Float64}
        r = q[1:3]
        R = q[4:6]
        ddr = ustrip(u"m^3/s^2",G*M*1u"Msun")*(R-r)/max(norm(R-r),rₘ)^3
        ddR = [0,0,0]
        append!(ddr,ddR)
        return ddr
    end
    Integrator=[]
    for i ∈ 1:Nₛ
        append!(vₛ[i],V₀)
        append!(rₛ[i],R₀)
    end
    for i ∈ 1:Nₛ
        push!(Integrator,init(SecondOrderODEProblem(Gravity,vₛ[i],rₛ[i],(0,dT*N_steps)),dt=dT,McAte2()))
    end 
    snapshot = zeros(frames,Nₛ+1,3)
    time = zeros(frames)
        for i ∈ 1:frames
            for j ∈ 1:Nₛ
                snapshot[i,j,:] =Integrator[j].u[7:9]*3.240779289444365e-20
            end
            snapshot[i,Nₛ+1,:] = Integrator[1].u[10:12]*3.240779289444365e-20
            time[i]=Integrator[1].t

            for j ∈ 1:divs
                step!.(Integrator)
            end
        end
    return snapshot,time 
end

function gif_galaxy(data,time,Nₛ,xlimit=[0,0],ylimit=[0,0],zlimit=[0,0])
    if ((norm(xlimit) == 0.0) || (norm(ylimit) == 0.0) || (norm(zlimit) == 0.0))
        xlimit = [minimum(data[:,:,1]),maximum(data[:,:,1])]
        ylimit = [minimum(data[:,:,2]),maximum(data[:,:,2])]
        zlimit = [minimum(data[:,:,3]),maximum(data[:,:,3])]
    end
    @gif for i ∈ 1:length(time)
        scatter3d(xlim=xlimit,ylim=ylimit,zlim=zlimit)
        for j ∈ 1:Nₛ
            scatter3d!(data[i:i,j,1],data[i:i,j,2],data[i:i,j,3], color=:red,legends=false,markersize=1.5)
           # i > trail ? plot3d!(Rotated[j][1,1:i],Rotated[j][2,1:i],Rotated[j][3,1:i]) : nothing
        end
        scatter3d!(data[i:i,Nₛ+1,1],data[i:i,Nₛ+1,2],data[i:i,Nₛ+1,3], color=:black,legends=false)
    end 
end

function gif_two_galaxies(data,time,N₁,N₂,xlimit=[0,0],ylimit=[0,0],zlimit=[0,0])
    if ((norm(xlimit) == 0.0) || (norm(ylimit) == 0.0) || (norm(zlimit) == 0.0))
        xlimit = [minimum(data[1:40,:,1]),maximum(data[1:40,:,1])]
        ylimit = [minimum(data[1:40,:,2]),maximum(data[1:40,:,2])]
        zlimit = [minimum(data[1:40,:,3]),maximum(data[1:40,:,3])]
    end
    @gif for i ∈ 1:length(time)
        scatter3d(xlim=xlimit,ylim=ylimit,zlim=zlimit)
        for j ∈ 1:N₁
            scatter3d!(data[i:i,j,1],data[i:i,j,2],data[i:i,j,3], color=:red,legends=false,markersize=1.5)
        end
        for j ∈ 1:N₂
            scatter3d!(data[i:i,N₁+j,1],data[i:i,N₁+j,2],data[i:i,N₁+j,3], color=:blue,legends=false,markersize=1.5)
           # i > trail ? plot3d!(Rotated[j][1,1:i],Rotated[j][2,1:i],Rotated[j][3,1:i]) : nothing
        end
        scatter3d!(data[i:i,N₁+N₂+1,1],data[i:i,N₁+N₂+1,2],data[i:i,N₁+N₂+1,3], color=:black,legends=false)
        scatter3d!(data[i:i,N₁+N₂+2,1],data[i:i,N₁+N₂+2,2],data[i:i,N₁+N₂+2,3], color=:black,legends=false)
    end 
end

function evolve_two_disks(primary,secondary,dT=1e-4u"yr",N_steps=100000,frames=500)
    #Integration stuff
    divs = trunc(Int64,N_steps/frames)
    dT=ustrip(u"s",dT)
    r₁ₘ,r₂ₘ=ustrip(u"m",primary["softening"]*primary["radius"]),ustrip(u"m",secondary["softening"]*secondary["radius"])
    N₁ₛ,N₂ₛ=primary["Nₛ"],secondary["Nₛ"]
    
    #Galaxy stuff
    M = [ustrip(u"Msun",primary["mass"]),ustrip(u"Msun",secondary["mass"])]
    R₁,R₂ = ustrip.(u"m",primary["center_pos"]),ustrip.(u"m",secondary["center_pos"])
    V₁,V₂ = ustrip.(u"m/s",primary["center_vel"]),ustrip.(u"m/s",secondary["center_vel"])

    #star stuff
    r₁,r₂=primary["star_pos"],secondary["star_pos"]
    v₁,v₂=primary["star_vel"], secondary["star_vel"]


    function Gravitus(dq::Vector{Float64},q::Vector{Float64},m::Vector{Float64},t::Float64)::Vector{Float64}
        pegnor= q[1:3]
        R₁,R₂= q[4:6],q[7:9]
        ddpegnor = ustrip(u"m^3/s^2",G*m[1]*1u"Msun")*(R₁-pegnor)/max(norm(R₁-pegnor),r₁ₘ)^3+ustrip(u"m^3/s^2",G*m[2]*1u"Msun")*(R₂-pegnor)/max(norm(R₂-pegnor),r₂ₘ)^3
        ddR₁ = ustrip(u"m^3/s^2",G*m[1]*1u"Msun")*(R₂-R₁)/max(norm(R₁-R₂),r₁ₘ)^3
        ddR₂ = -ddR₁
        append!(ddpegnor,ddR₁,ddR₂)
        return ddpegnor
    end
    
    Integrator=[]
    for i ∈ 1:N₁ₛ
        append!(v₁[i],V₁,V₂)
        append!(r₁[i],R₁,R₂)
    end
    for i ∈ 1:N₂ₛ
        append!(v₂[i],V₁,V₂)
        append!(r₂[i],R₁,R₂)
    end
    
    for k ∈ 1:N₁ₛ
        push!(Integrator,init(SecondOrderODEProblem(Gravitus,v₁[k],r₁[k],(0,N_steps*dT),M),dt=dT,McAte2()))
    end
    
    for i ∈ 1:N₂ₛ
        push!(Integrator,init(SecondOrderODEProblem(Gravitus,v₂[i],r₂[i],(0,N_steps*dT),M),dt=dT,McAte2()))
    end 
    snapshot = zeros(frames,N₁ₛ+N₂ₛ+2,3)
    time = zeros(frames)
        for i ∈ 1:frames
            for j ∈ 1:N₁ₛ
                snapshot[i,j,:] =Integrator[j].u[10:12]*3.240779289444365e-20
            end
            for j ∈ 1:N₂ₛ
                snapshot[i,N₁ₛ+j,:] =Integrator[N₁ₛ+j].u[10:12]*3.240779289444365e-20
            end
            snapshot[i,N₁ₛ+N₂ₛ+1,:] = Integrator[1].u[13:15]*3.240779289444365e-20
            snapshot[i,N₁ₛ+N₂ₛ+2,:] = Integrator[1].u[16:18]*3.240779289444365e-20
            time[i]=Integrator[1].t * ustrip(u"Myr",1u"s")
            for j ∈ 1:divs
                step!.(Integrator)
            end
        end
    return snapshot,time 
end