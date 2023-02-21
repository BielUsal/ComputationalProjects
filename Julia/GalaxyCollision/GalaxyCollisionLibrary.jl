using LinearAlgebra
using Unitful
using UnitfulAstro

G = 4.3009E-3 *1u"pc *(km/s)^2 /Msun"
function format_parameters(galaxy_args)
    #I didn't have this at first, but this makes all galaxies have an uniform format. Not too crazy
    return Dict(
        "mass" => galaxy_args[1]*1u"Msun",
        "radius"     => galaxy_args[2]*1u"pc",
        "center_pos" => galaxy_args[3].*1u"pc",
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
        u = normalize(u)

        Rotation = [
        u[1]*u[1]*(1-cosθ)+cosθ u[1]*u[2]*(1-cosθ)-u[2]*sinθ u[1]*u[3]*(1-cosθ)+u[1]*sinθ;

        u[2]*u[1]*(1-cosθ)+u[3]*sinθ u[2]*u[2]*(1-cosθ)+cosθ u[2]*u[3]*(1-cosθ)-u[1]*sinθ;

        u[3]*u[1]*(1-cosθ)+u[2]*sinθ u[3]*u[1]*(1-cosθ)+u[1]*sinθ u[3]*u[3]*(1-cosθ)+cosθ
        ]
    end
#-----------------------------------------------------------------------------------------------------------------------------------
    galaxy["stars_pos"] = []
    galaxy["stars_vel"] = []

    R = galaxy["softening"] * galaxy["radius"]
    for i ∈ 1:galaxy["Nᵣ"] #iterate over the rings
        #Randomly distribute the initial radii and angles 
        rₛ =(R * ones(StarsPerRing) + dr *rand(StarsPerRing))./1u"pc" #creates an n-dimensional vector of radii from R to R+dr 
        ϕₛ = 2π *rand(StarsPerRing)

        # Positions 
        vecᵣ = (Rotation * ([rₛ.*cos.(ϕₛ),rₛ.*sin.(ϕₛ),zeros(StarsPerRing)])).*1u"pc"
        x = galaxy["center_pos"][1].+vecᵣ[1]
        y = galaxy["center_pos"][2].+vecᵣ[2]
        z = galaxy["center_pos"][3].+vecᵣ[3]

        # Velocities
        Tₛ = 2π * uconvert.(u"s",sqrt.((rₛ*1u"pc").^3/(G *Introoder["mass"])))

        Δϕ = 2π *uconvert(u"s",dT)./Tₛ 

        vecᵥ = (Rotation* [(rₛ/(uconvert(u"s",dT)*1u"1/s")).*(cos.(ϕₛ)-cos.(ϕₛ-Δϕ)),rₛ/(uconvert(u"s",dT)*1u"1/s").*(sin.(ϕₛ)-sin.(ϕₛ-Δϕ)),zeros(StarsPerRing)])*1u"pc/s"
        v₁ = uconvert.(u"km/s",galaxy["center_vel"][1].+vecᵥ[1])
        v₂ = uconvert.(u"km/s",galaxy["center_vel"][2].+vecᵥ[2])
        v₃ = uconvert.(u"km/s",galaxy["center_vel"][3].+vecᵥ[3])
        push!(galaxy["stars_pos"],[x,y,z])
        push!(galaxy["stars_vel"],[v₁,v₂,v₃])
        R += dr
    end 
    galaxy["vel_scale"]=uconvert(u"km/s",√(G*galaxy["mass"]/(0.5*R)))
end