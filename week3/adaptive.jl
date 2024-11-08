using Plots

function doublepend(x,p,t)
    # parameters
    m1,m2,l1,l2, g= p

    dx = zeros(length(x))

    ϕ1 = x[1]
    ϕ2 = x[2]
    q1 = x[3]
    q2 = x[4]

    cosphidiff = cos(ϕ1-ϕ2)
    sinphidiff = sin(ϕ1-ϕ2)
    detM = l1^2*l2^2*m2*(m1+m2) - (m2*l1*l2*cosphidiff)^2
    
    # ODE
    dx[1] = (q1*m2*l2^2 - q2*m2*l1*l2*cosphidiff)/detM
    dx[2] = (-q1*m2*l1*l2*cosphidiff + q2*l1^2*(m1+m2))/detM
    #dx[2] = ( q2*(m1+m2)/(m2*l2^2) - q1/(l1*l2)*cosphidiff ) * ( m1+m2-cosphidiff^2*m2 )^(-1)
    #dx[1] = q1/((m1+m2)*l1^2) - m2/(m1+m2) * l2/l1 * dx[2]*cosphidiff
    dx[3] = - m2*l1*l2*dx[1]*dx[2]*sinphidiff - g*l1*sin(ϕ1)*(m1+m2)
    dx[4] = m2*l1*l2*dx[1]*dx[2]*sinphidiff - g*l2*sin(ϕ2)*m2

    return dx#,ϕ1,ϕ2
end

function rad(deg)
    return deg*pi/180
end

function to_xy(phi1,phi2)
    x1 = l1*sin.(phi1)
    y1 = -l1*cos.(phi1)
    x2 = x1 + l2*sin.(phi2)
    y2 = y1 - l2*cos.(phi2)
    return (x1,y1,x2,y2)
end

function rk4(x,p,t,h)
    k1 = doublepend(x,p,t)
    k2 = doublepend(x+0.5*h*k1,p,t+0.5*h)
    k3 = doublepend(x+0.5*h*k2,p,t+0.5*h)
    k4 = doublepend(x+h*k3,p,t+h)
    x = x + 1/6*h*(k1 + 2*k2 + 2*k3 + k4)
    return x#1/6*(k1 + 2*k2 + 2*k3 + k4)#x_out
end


function energy(phi1,phi2,q1,q2,p)
    m1,m2,l1,l2, g = p

    cosphidiff = cos(phi1 - phi2)
    sinphidiff = sin(phi1 - phi2)

    phi2dot = (q2 * (m1 + m2) / (m2 * l2^2) - q1 / (l1 * l2) * cosphidiff) / (m1 + m2 - cosphidiff^2 * m2)
    phi1dot = q1 / ((m1 + m2) * l1^2) - m2 / (m1 + m2) * (l2 / l1) * phi2dot * cosphidiff

    kinetic = 0.5*m1*(l1*phi1dot)^2 + 0.5*m2*( (l1*phi1dot)^2+(l2*phi2dot)^2 + 2*l1*l2*phi1dot*phi2dot*cos(phi1-phi2) )
    potential = m1*g*l1*(1-cos(phi1)) + m2*g*(l1*(1-cos(phi1))+l2*(1-cos(phi2)))
    #potential = -m1*g*l1*cos(phi1) - m2*g*(l1*cos(phi1)+l2*cos(phi2))
    E = kinetic + potential
    return E
end


x0 = [rad(50), rad(-120),0,0]

m1,m2,l1,l2,g = 0.5, 1.0,2.0,1.0,1.0

params = [m1,m2,l1,l2,g]

h = 0.01 # 1e-4
h_max = h
t_max = 100

# allocate vectors
phi1_arr = Float64[]#zeros(timesteps)
phi2_arr = Float64[]#zeros(timesteps)
q1_arr = Float64[]#zeros(timesteps)
q2_arr = Float64[]#zeros(timesteps)
E = Float64[]
t_arr =  Float64[]
h_vals = Float64[]

push!(h_vals,h)


# initializing
push!(phi1_arr,x0[1])
push!(phi2_arr,x0[2])
push!(q1_arr,x0[3])
push!(q2_arr,x0[4])
E0 = energy(phi1_arr[1],phi2_arr[1],0,0,params)
println(E0)
push!(E,E0)

epsilon = 1e-8

reduce_factor = 0.5
increase_factor = 1.2
h_min = 1e-4

i = 2
t_current = 0
push!(t_arr,t_current)
while t_current <= t_max
    # solving ODE
    old_state = [phi1_arr[i-1],phi2_arr[i-1],q1_arr[i-1],q2_arr[i-1]]
    new_state = rk4(old_state,params,t_arr[i-1],h)

    #phi1_arr[i],phi2_arr[i],q1_arr[i],q2_arr[i] = new_state
    push!(phi1_arr,new_state[1])
    push!(phi2_arr,new_state[2])
    push!(q1_arr,new_state[3])
    push!(q2_arr,new_state[4])

    # energy and step size control
    E_t = energy(phi1_arr[i],phi2_arr[i],q1_arr[i],q2_arr[i],params)
    if abs(E_t-E[i-1])/E[i-1] > epsilon
        global h *= reduce_factor
    elseif abs(E_t-E[i-1])/E[i-1] < epsilon/10 # conservative check for when to increase_factor
        global h *= increase_factor
    end
    if h < h_min
        global h = h_min
    end
    if h > h_max
        global h = h_max
    end
    push!(h_vals,h)
    push!(E,E_t)
    global t_current += h
    push!(t_arr,t_current)


    if i%1000 == 0
        println(t_current)
    end
    global i += 1
end

timesteps = length(t_arr)
x1,y1,x2,y2 = zeros(timesteps),zeros(timesteps),zeros(timesteps),zeros(timesteps)
x1,y1,x2,y2 = to_xy(phi1_arr,phi2_arr)

p_pos = plot(x1,y1,xlabel="x",ylabel="y")
plot!(p_pos,x2,y2)

p_angles = plot(t_arr,phi1_arr,xlabel="time",ylabel="angle")
plot!(p_angles,t_arr,phi2_arr)

E = abs.(E.-E0)./E0

pE = plot(t_arr,E,xlabel="time",ylabel="ΔE")

ph = plot(t_arr,h_vals,xlabel="time",ylabel="timestep")

fig = plot(
    p_pos, p_angles, pE, ph,
    layout=(2, 2),  # 2x2 grid layout
    size=(800, 600) # Adjust the figure size if needed
)


#savefig(fig,"adaptive$epsilon+red.png")
#savefig(fig,"adaptive$epsilon+lowstart.png")
savefig(fig,"adaptive$epsilon+prev_lowstart.png")
#savefig(fig,"plots_adaptive/not_adaptive_$h+_h.png")