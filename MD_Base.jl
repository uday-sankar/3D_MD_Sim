# Base Functions
# Velocity-verlet algorithm is used for time propogation
# Temeprature is defiend as average Kinetic enenrgy
# All of this code assumes we particles in 3D space.
# But many of the functions defined here work for N-degrees of freedom.
using Plots
using LinearAlgebra

Har(xyz,eq=2*ones(3)) = 1/2*sum((xyz - eq).^2)
F_Har(x,y) = [2-x,2-y]

##Damped trajectory
function Run_damp(XYZ,Vxyz;T=1:0.1:10,F=F_Har,V=Har,m=-1,bl=20)
    #INPUT: initial cordinates and velocities and the time axis and number of particles
    #RETURNS: The X,Y as a function of time. Total energy as a function of time
    #damping works by scaling down (here making it zero) the velocity whenever aceleration becomes negative (i.e a barrier)
    print("\nRunning a Damped Trajectory")
    xyz_dim = size(XYZ)
    #Sanity check
    if xyz_dim != size(Vxyz)
        print("!!Dimension Mismatch: XYZ and Velocity matrix has different dimesnions")
        exit()
    elseif m==-1
        m = ones(xyz_dim[1])
    end 
    np = xyz_dim[1] #length(Xinit)
    dt = T[2] - T[1] # Assuming uniform grid
    lent = length(T)
    XYZ_traj = zeros(lent, xyz_dim...) #The XYZ trajectory
    E, Temp, PE = zeros(lent), zeros(lent), []
    XYZ_new, Vxyz_new = copy(XYZ), copy(Vxyz)
    Ft = zeros(xyz_dim)
    PE_int = 0
    for i in 1:np
        for j in i+1:np
            ft = F(XYZ_new[i,:]-XYZ_new[j,:])
            PE_int += V(XYZ_new[i,:]-XYZ_new[j,:])
            Ft[i, :] += ft
            Ft[j, :] -= ft
        end
    end
    KE = sum( m'*(Vxyz_new.^2) )/2.0
    append!(PE, PE_int)
    E[1] = KE + PE_int
    Temp[1] = KE/np
    #print("\nStep 1 PE = $(PE[1])")
    for i in 2:lent
        XYZ_new = XYZ_new + Vxyz_new*dt + 0.5*Ft./m*dt^2 
        Ftn = zeros(xyz_dim)
        PE_int = 0 #Potential energy due to inter-particle interactions
        for i in 1:np-1
            for j in i+1:np
                ftn = F(XYZ_new[i,:]-XYZ_new[j,:])
                PE_int += V(XYZ_new[i,:]-XYZ_new[j,:])
                Ftn[i, :] += ftn/m[i]
                Ftn[j, :] -= ftn/m[j]
            end
        end
        F_mid = (Ftn + Ft)/2.0 
        Vxyz_new = Vxyz_new + F_mid*dt
        KE = sum(m'*(Vxyz_new.^2))/2.0
        XYZ_traj[i, :, :] = XYZ_new
        append!(PE,PE_int)
        E[i], Temp[i]   =  KE + PE_int,  KE/np
        Ft = copy(Ftn)
        for j in 1:np
            v_vec = Vxyz_new[j,:]
            F_vec = F_mid[j,:]
            F_vec = F_vec./norm(F_vec)
            v_proj = v_vec'F_vec
            if v_proj  < 0
                Vxyz_new[j,:] = zeros(xyz_dim[2])
            else
                Vxyz_new[j,:] = v_proj.*F_vec
            end
        end
        print("\n\tStep $i PE = $(PE[i])")
        if ((abs(PE[i]-PE[i-1]) < 1e-15) && (sum(abs.(Ft)) < 1e-12) && sum(abs.(Vxyz_new)) < 1e-12 )
            print("\nConvergence reached after: $i steps")
            break
        end
    end
    return XYZ_traj, E, Temp, PE
end

function Run(XYZ,Vxyz;T=1:0.1:10,F=F_Har,V=Har,m=-1,bl=20)
    #INPUT: initial cordinates and velocities and the time axis and number of particles
    #RETURNS: The X,Y as a function of time. Total energy as a function of time
    print("\n Running Simulation")
    xyz_dim = size(XYZ)
    #Sanity check
    if xyz_dim != size(Vxyz)
        print("!!Dimension Mismatch: XYZ and Velocity matrix has different dimesnions")
        exit()
    elseif m == -1
        m = ones(xyz_dim[1])
    end 
    np = xyz_dim[1] #length(Xinit)
    print("\n \tNumber of particles:$np")
    dt = T[2] - T[1] # Assuming uniform grid
    lent = length(T)
    XYZ_traj = zeros(lent, xyz_dim...) #The XYZ trajectory
    E, Temp, PE = zeros(lent), zeros(lent), zeros(lent)
    XYZ_new, Vxyz_new = copy(XYZ), copy(Vxyz)
    Ft = zeros(xyz_dim)
    PE_int = 0
    for i in 1:np
        for j in i+1:np
            ft = F(XYZ_new[i,:]-XYZ_new[j,:])
            PE_int += V(XYZ_new[i,:]-XYZ_new[j,:])
            Ft[i, :] += ft/m[i]
            Ft[j, :] -= ft/m[j]
        end
    end
    KE = sum(m'*(Vxyz_new.^2))/2.0
    PE[1] = PE_int
    E[1] = KE + PE_int
    Temp[1] = KE/np
    for i in 2:lent
        XYZ_new = XYZ_new + Vxyz_new*dt + 0.5*Ft*dt^2 
        Ftn = zeros(xyz_dim)
        PE_int = 0 #Potential energy due to inter-particle interactions
        for i in 1:np-1
            for j in i+1:np
                ftn = F(XYZ_new[i,:]-XYZ_new[j,:])
                PE_int += V(XYZ_new[i,:]-XYZ_new[j,:])
                Ftn[i, :] += ftn/m[i]
                Ftn[j, :] -= ftn/m[j]
            end
        end
        F_mid = (Ftn + Ft)/2.0 
        Vxyz_new = Vxyz_new + F_mid*dt
        KE = sum(m'*(Vxyz_new.^2))/2.0
        XYZ_traj[i, :, :] = XYZ_new
        PE[i], E[i], Temp[i]   = PE_int, KE + PE_int,  KE/np
        cond = abs.(XYZ_new) .>= bl #Condition for bounce from boundary
        Vxyz_new[cond] .= -Vxyz_new[cond] #Velocity inversion for bounce
        Ft, Ft = copy(Ftn), copy(Ftn)
    end
    return XYZ_traj, E, Temp, PE
end

function Grid_Initialize(N=[5,5,5],d=[2,2,2])
    # Ng is the number of grid points in x and y direction
    # d is the distance between two images
    Nx, Ny, Nz = N
    dx, dy, dz = d 
    ##
    set_x =  -(Nx-1)*dx/2:dx:(Nx-1)*dx/2
    set_y =  -(Ny-1)*dy/2:dy:(Ny-1)*dy/2
    set_z =  -(Nz-1)*dz/2:dz:(Nz-1)*dz/2
    np = length(set_x)*length(set_y)*length(set_z)
    XYZ = zeros(np,3)
    #particle index
    n=1
    for i in 1:Nx
        for j in 1:Ny
            for k in 1:Nz
                XYZ[n,:] = [set_x[i], set_y[j], set_z[k]]
                n+=1
            end
        end
    end
    ##
    Vxyz = zeros(np,3)
    return XYZ, Vxyz
end

function G_Stab(Np,dp;T = 0:0.001:400,F=Har,V=F_Har,m=-1,sim=0) #Short form for Grid stabilization
    #Np is number of particles across each grid axis
    print("Setting Initial Conditions")
    dt = T[2] - T[1]
    #lt = length(T)
    ##  Grid Initialization of position coordinates 
    print("\nSetting initial conditions for stabilizing")
    XYZ, Vxyz = Grid_Initialize(Np,dp)
    num_particles = size(XYZ)[1]
    if m == -1
        m = ones(size(XYZ)[1])
    end
    ##
    print("\nNo of Particles:$num_particles")
    scatter(XYZ[:,1],XYZ[:,2],XYZ[:,3],legend = false,color=:blue,linewidth=1)
    savefig("GStab_Init.png")
    ##
    print("\nRunning Stabilization") 
    XYZt, Et, Te, PE = Run_damp(XYZ,Vxyz;T=T,V=V,F=F,m=m)
    print("\nChange in Potential Energy:",PE[end]-PE[1])
    plot(T[1:length(PE)],PE,xlabel="T",ylabel="V")
    savefig("Gstab_V.png")
    ##
    if sim == 1
        Traj_animation_3D(XYZt[1:length(PE),:,:];filename="Gstab.gif",frames=200)
    end
    ##
    XYZend = copy(XYZt[length(PE),:,:]) 
    scatter(XYZend[:,1],XYZend[:,2],XYZend[:,3],legend = false,color=:blue,linewidth=1)
    savefig("Gstab_Stable.png")
    plot(T[1:length(PE)],Te[1:length(PE)])
    savefig("Gstab_Temp.png")
    CoM = (m'*XYZend)/num_particles
    XYZ_norm = XYZend .- CoM
    open("Stabilized_geom.dat","w") do f
        for i in 1:num_particles
            write(f,"$i $(XYZ_norm[i,1]) $(XYZ_norm[i,2]) $(XYZ_norm[i,3])\n" )
        end
    end
    return XYZ_norm
end

function R_Stab(Np,xyz_range;T = 0:0.001:400,F=Har,V=F_Har,m=-1,sim=0,fold=".") #Short form for Random stabilization
    # Np  is the total number of particles
    # xyz_range is the range in which the the random numbers are to be initialized
    # It should be an vector of size 3
    #dt = T[2] - T[1]
    #lt = length(T)
    ##  Grid Initialization of position coordinates 
    print("\n Random Stabilization Started")
    print("\n \tNo of Particles:$num_particles")
    print("\n \tSetting initial conditions for random stabilizing")
    XYZ = hcat(xyz_range[1]*rand(Np),xyz_range[2]*rand(Np),xyz_range[3]*rand(Np))#rand(Np,3)
    Vxyz = zeros(Np,3)
    if m == -1
        m = ones(size(XYZ)[1])
    end
    ##
    scatter(XYZ[:,1],XYZ[:,2],XYZ[:,3],legend = false,color=:blue,linewidth=1)
    savefig("$fold/RStab_Init.png")
    ##
    print("\n \tRunning Stabilization") 
    XYZt, Et, Te, PE = Run_damp(XYZ,Vxyz;T=T,V=V,F=F,m=m)
    print("\n \tChange in Potential Energy:",PE[end]-PE[1])
    plot(T[1:length(PE)],PE,xlabel="T",ylabel="V")
    savefig("$fold/RStab_V.png")
    ##
    if sim == 1
        Traj_animation_3D_stab(XYZt[1:length(PE),:,:];filename="$fold/Rstab.gif",frames=200)
    end
    ##
    XYZend = copy(XYZt[length(PE),:,:]) 
    scatter(XYZend[:,1],XYZend[:,2],XYZend[:,3],legend = false,color=:blue,linewidth=1)
    savefig("$fold/RStab_Stable.png")
    plot(T[1:length(PE)],Te[1:length(PE)])
    savefig("$fold/RStab_Temp.png")
    CoM = (m'*XYZend)/num_particles
    XYZ_norm = XYZend .- CoM
    open("$fold/RStab_geom.dat","w") do f
        for i in 1:num_particles
            write(f,"$i $(XYZ_norm[i,1]) $(XYZ_norm[i,2]) $(XYZ_norm[i,3])\n" )
        end
    end
    print("\n \tRandom Stabilization Done")
    return XYZ_norm
end

function G_Run(Np,dp;F=Har,V=F_Har) #Short form for Grid stabilization
    print("Setting Initial Conditions")
    dt =0.0001
    T = 0:dt:1
    #lt = length(T)
    m = ones(np)
    ##  Grid Initialization of position coordinates 
    Xinit, Yinit, Vxinit, Vyinit = Grid_Initialize(Np,dp)
    ##
    scatter(Xinit,Yinit,legend=false,color=:blue,linewidth=1)
    savefig("GRun_Init.png")
    ##
    print("\nRunning Simulation") 
    Xt,Yt,Et, Te, PE = Run(Xinit,Yinit,Vxinit,Vyinit;T=T,V=V,F=F,m=m)
    print("\nEnergy Error:",maximum(Et)-minimum(Et))
    ##
    Xend, Yend = copy(Xt[end,:]), copy(Yt[end,:]) 
    scatter(Xend,Yend,legend=false,color=:blue,linewidth=1)
    savefig("GRun_Final.png")
    plot(T,Te)
    savefig("GRun_Temp.png")
    Xend = Xend .- (sum(Xend)/25)
    Yend = Yend .- (sum(Yend)/25)
    open("GRun_Final.dat","w") do f
        for i in 1:np
            write(f,"$i $(Xend[i]) $(Yend[i])\n" )
        end
    end
    return Xt, Yt
end

#function Traj_animation_3D(Traj;m=-1, frames=100, bl=20, filename="Test.gif")
#    if m == -1
#        m = ones(size(Traj)[2])
#    end 
#    print("\nAnimating\n")
#    lt = size(Traj)[1]
#    N_step_frame =  Int(ceil(lt/frames))
#    plot()
#    anim = @animate for j in 1:frames
#        i = Int(ceil(j*lt/frames))
#        plot([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], -bl.*ones(5), color = :black, linewidth=3, legend=false, xlims=(-bl, bl), ylims=(-bl, bl), zlims=(-bl, bl), xlabel="X", ylabel="Y", zlabel="Z")
#        plot!([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], bl.*ones(5), color = :black, linewidth=3, legend=false)
#        plot!(bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], color = :black, linewidth=3, legend=false)
#        plot!(-bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], color = :black, linewidth=3, legend=false)
#        #plot()
#        for j in 1:size(Traj)[2]
#            scatter!([Traj[i, j, 1]], [Traj[i, j, 2]], [Traj[i, j, 3]], color=:blue, markersize=2*m[j])
#        end
#        if (i-8*N_step_frame) > 0
#            plot!([Traj[(i-8*N_step_frame):i,: , 1]], [Traj[(i-8*N_step_frame):i,: , 2]], [Traj[(i-8*N_step_frame):i,: , 3]], color=:blue, linewidth=1)
#        #else
#        #    plot!([Traj[1:i,:,1]], [Traj[1:i,:,2]], [Traj[1:i,:,2]], color=:blue, linewidth=1)
#        end
#    end
#    gif(anim, filename, fps=30)
#end

function Traj_animation_3D(Traj;m=-1, frames=100, bl=20, filename="Test.gif")
    if m == -1
        m = ones(size(Traj)[2])
    end 
    print("\nAnimating\n")
    lt = size(Traj)[1]
    N_step_frame =  Int(ceil(lt/frames))
    anim = @animate for j in 1:frames
        i = Int(ceil(j*lt/frames))
        plot([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], -bl.*ones(5), color = :black, linewidth=3, legend=false, xlims=(-bl, bl), ylims=(-bl, bl), zlims=(-bl, bl), xlabel="X", ylabel="Y", zlabel="Z")
        plot!([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], bl.*ones(5), color = :black, linewidth=3, legend=false)
        plot!(bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], color = :black, linewidth=3, legend=false)
        plot!(-bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], color = :black, linewidth=3, legend=false)
        for j in 1:size(Traj)[2]
            scatter!([Traj[i, j, 1]], [Traj[i, j, 2]], [Traj[i, j, 3]], color=:blue, markersize=2*m[j])
        end
        if (i-8*N_step_frame) > 0
            plot!([Traj[(i-8*N_step_frame):i,: , 1]], [Traj[(i-8*N_step_frame):i,: , 2]], [Traj[(i-8*N_step_frame):i,: , 3]], color=:blue, linewidth=1)
        #else
        #    plot!([Traj[1:i,:,1]], [Traj[1:i,:,2]], [Traj[1:i,:,2]], color=:blue, linewidth=1)
        end
    end
    gif(anim, filename, fps=30)
end

function Traj_animation_3D_stab(Traj;m=-1, frames=100, bl=20, filename="Test.gif")
    if m == -1
        m = ones(size(Traj)[2])
    end 
    print("\nAnimating\n")
    lt = size(Traj)[1]
    N_step_frame =  Int(ceil(lt/frames))
    anim = @animate for j in 1:frames
        i = Int(ceil(j*lt/frames))
        for j in 1:size(Traj)[2]
            scatter([Traj[i, j, 1]], [Traj[i, j, 2]], [Traj[i, j, 3]], color=:blue, xlabel="X", ylabel="Y", zlabel="Z",legend=false)
        end
        if (i-8*N_step_frame) > 0
            plot!([Traj[(i-8*N_step_frame):i,: , 1]], [Traj[(i-8*N_step_frame):i,: , 2]], [Traj[(i-8*N_step_frame):i,: , 3]], color=:blue, linewidth=1)
        #else
        #    plot!([Traj[1:i,:,1]], [Traj[1:i,:,2]], [Traj[1:i,:,2]], color=:blue, linewidth=1)
        end
    end
    gif(anim, filename, fps=30)
end

function Traj_Plot_new(Traj; frames=100, bl=20, filename="Test.gif")
    println("\nAnimating\n")
    plot()
    lt = size(Traj)[1]
    N_step_frame =  Int(ceil(lt/frames))
    anim = @animate for j in 1:frames
        i = Int(ceil(j*lt/frames))
        #plot!([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], zeros(5), color = :black, linewidth=3, legend=false, xlims=(-25, 25), ylims=(-25, 25), zlims=(-25, 25), xlabel="X", ylabel="Y", zlabel="Z")
        plot( scatter(
            x=[Traj[i,:,1]],
            y=[Traj[i,:,2]],
            z=[Traj[i,:, 3]],
            mode="markers",
            marker=attr(
                size=12,
                color=:blue,                # set color to an array/list of desired values
                #colorscale="Viridis",   # choose a colorscale
                opacity=0.8
            ),
            type="scatter3d"
        ), 
        Layout(margin=attr(l=0, r=0, b=0, t=0)))
        #if j > 10
        #    plot!([Traj[(i-4*N_step_frame):i,: , 1]], [Traj[(i-4*N_step_frame):i,: , 2]], [Traj[(i-4*N_step_frame):i,: , 3]], color=:blue, linewidth=1)
        #else
        #    #scatter!([Traj[i,:,1]], [Traj[i,:,2]], [Traj[i,:,3]], color=:blue, markersize=2)
        #    plot!([Traj[1:i,:,1]], [Traj[1:i,:,2]], [Traj[1:i,:,2]], color=:blue, linewidth=1)
        #end
    end
    gif(anim, filename, fps=30)
end

function Plot_frame(XYZ; file_name="Test_frame.png",m=-1,bl=20)
    print("\nPlotting frame")
    if m == -1
        m = ones(size(XYZ)[1])
    end 
    plot([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], -bl.*ones(5), color = :black, linewidth=3, legend=false, xlims=(-bl, bl), ylims=(-bl, bl), zlims=(-bl, bl), xlabel="X", ylabel="Y", zlabel="Z")
    plot!([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], bl.*ones(5), color = :black, linewidth=3, legend=false)
    plot!(bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], color = :black, linewidth=3, legend=false)
    plot!(-bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], color = :black, linewidth=3, legend=false)
    #plot()
    for i in 1:size(XYZ)[1]
        scatter!([XYZ[i,1]], [XYZ[i,2]], [XYZ[i, 3]], color=:blue, markersize=2*m[i])
    end
    savefig(file_name)
end

#Change in Potential Energy:-66.82651854846782 With break conditions of 1e-15