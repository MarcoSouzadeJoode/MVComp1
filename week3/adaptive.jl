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
    dx[3] = - m2*l1*l2*dx[1]*dx[2]*sinphidiff - g*l1*sin(ϕ1)*(m1+m2)
    dx[4] = m2*l1*l2*dx[1]*dx[2]*sinphidiff - g*l2*sin(ϕ2)*(m1+m2)

    return dx,x
end

function rk4(x,p,t,h)
    k1,phis = doublepend(x,p,t)
    k2,phis = doublepend(x+0.5*h*k1,p,t+0.5*h)
    k3,phis = doublepend(x+0.5*h*k2,p,t+0.5*h)
    k4,phis = doublepend(x+h*k2,p,t+h)
    x = x + 1/6*h*(k1 + 2*k2 + 2*k3 + k4)
    return x
end

function energy(x1,y1,x2,y2,p,timerange)
    m1,m2,l1,l2, g= p
    E = zeros(length(timerange))
    E0 = 0
    for t in 1:length(timerange)-1
        kinetic = 1/2*( m1*(((x1[t]-x1[t+1])/h)^2 + ((y1[t]-y1[t+1])/h)^2) + m2*(((x2[t]-x2[t+1])/h)^2 + ((y2[t]-y2[t+1])/h)^2))
        potential = - m1*g*y1[t] - m2 * g*y2[t]
        E[t] = kinetic + potential
        if t == 2 
            E0 = E[t]
        end
    end 
    return E, E0
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

x0 = [rad(50), rad(-120),0,0]

m1,m2,l1,l2,g = 0.5, 1.0,2.0,1.0,1.0

params = [m1,m2,l1,l2,g]

h = 0.05

t_range = 0:h:100

timesteps = length(t_range)

# allocate vectors
phi1_arr = zeros(timesteps)
phi2_arr = zeros(timesteps)
q1_arr = zeros(timesteps)
q2_arr = zeros(timesteps)

phi1_arr[1],phi2_arr[1],q1_arr[1],q2_arr[1] = x0

for i in 2:timesteps
    old_state = [phi1_arr[i-1],phi2_arr[i-1],q1_arr[i-1],q2_arr[i-1]]
    new_state = rk4(old_state,params,t_range[i],h)
    phi1_arr[i],phi2_arr[i],q1_arr[i],q2_arr[i] = new_state
end


x1,y1,x2,y2 = zeros(timesteps),zeros(timesteps),zeros(timesteps),zeros(timesteps)
x1,y1,x2,y2 = to_xy(phi1_arr,phi2_arr)

p = plot(x1,y1)
plot!(p,x2,y2)
display(p)

E,E0 = energy(x1,y1,x2,y2,params,t_range)

E = abs.(E.-E0)./E0

pE = plot(t_range,E)
display(pE)