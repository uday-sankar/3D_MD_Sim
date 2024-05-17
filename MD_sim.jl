#A set of 2D interacting particles. The interaction is modelled using a Morse Potential.
#The particles are initalized from auniform grid. 
#using Plots
#using PlotlyJS
#plotlyjs()
##
include("MD_Base.jl")
##
function LJ_v(xyz;ϵ=Energy_param,σ=Sig)# Defing the interaction between particles to be LJ potential
    #ϵ = epsilon parameter, σ is the sigma parametre defined
    r = (xyz'xyz)^0.5#root of (Distance between the particles squared)
    #defiing LJ potential
    V_lj =  4ϵ*((σ/r)^12 - (σ/r)^6)
    # A more computationally cheap approach
    #   r_6 = (xyz'xyz)^3#the internuclear dsitance raised to 6, r^6 
    #   V_lj = 4ϵ*(σ^6/r_6)*(σ^6/r_6 - 1)
    return V_lj
end
##
#Harmonic intearction is assumed between particles
function LJ_f(xyz;ϵ=Energy_param,σ=Sig)
    r = (xyz'xyz)^0.5#Distance between the particles 
    F = 24ϵ*(2*(σ/r)^12 - (σ/r)^6)/r^2
    F = F.*xyz
    return F
end
##
a_global = 1
De_global = 1
r_eq_global = 1.5
function V_Morse(xyz;a=a_global,De=De_global,r_eq=r_eq_global)# Defing the interaction between particles to be harmonic 
    r = (xyz'xyz)^0.5#root of (Distance between the particles squared)
    v = De*(1-exp(-a*(r-r_eq)))^2#Morse potential in 2 dimensions
    return v
end
##
#Harmonic intearction is assumed between particles
function F_Morse(xyz;a=a_global,De=De_global,r_eq=r_eq_global)
    r = (xyz'xyz)^0.5#root of (Distance between the particles squared)
    F = -2*De*(1-exp(-a*(r-r_eq)))*exp(-a*(r-r_eq))*a/r
    F = F.*xyz
    #Fx, Fy = F*xl/r, F*yl/r
    return F
end
##
##Interaction Parameters
Energy_param = 1
Sig = 3.0
##
a_global = 1
De_global = 1
r_eq_global = 3
#number of images in x direction 
X_num, dx = 3, 4.0 
Y_num, dy = 3, 4.0
Z_num, dz = 3, 4.0
Np = [X_num, Y_num, Z_num]
dp = [dx, dy, dz]
num_particles = X_num*Y_num*Z_num
##box parameters
bl = 15 # half-box length
##
F_global = LJ_f#F_Morse#LJ_f#F_Morse#LJ_xy
V_global = LJ_v#V_Morse#LJ_v#V_Morse#LJ_p
##
fold="Random"
##
T_stab = 0:0.05:400
#Xyz_stab = R_Stab(num_particles,[20,5,5];F=F_global,V=V_global,T=T_stab,sim=1,fold=fold)
Xyz_stab = G_Stab(Np,dp;F=F_global,V=V_global,T=T_stab,sim=1)#,fold=fold)
##
XYZ = Xyz_stab#20*rand(np,3).-10.0 # 4 partciles in 3 dimesnion
##
np = size(XYZ)[1]
m = ones(np)
#m[rand(1:np,10)] .= 5
##
V_XYZ = zeros(np,3)#0.0*rand(np,3)
V_XYZ[:,2] = ones(np)
##
Plot_frame(XYZ,file_name="GRun_Init.png")
#ß
T_run = 0:0.01:100
XYZ_traj,E,Temp,PE = Run(XYZ,V_XYZ;T=T_run,F=F_global,V=V_global,bl=bl,m=m)
##
p = plot(T_run,E,xlabel="T",ylabel="Total E")
savefig(p,"GRun_Total_energy.png")
print("\nEnergy Error:",maximum(E)-minimum(E))
p=plot(T_run,Temp)
savefig(p,"Run_Temp.png")
##
Traj_animation_3D(XYZ_traj,frames=400,filename="GRun.gif",bl=bl)
