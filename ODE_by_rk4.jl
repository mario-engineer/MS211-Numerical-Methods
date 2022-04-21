#=
Initial Value Problems
Let be the ODE: y'= y(x^2) - y; y(0)=1
Find an aproximate solution for 0 <= x <= 1
Use 4th Order Runge-Kutta Method with step h=0.1
The exact solution is y(x)=exp(-x + x^3/3)
=#
using Plots

t0=0.0
tf=1.0
h=0.1
n=convert(Int64, (tf-t0)/h)
y0=1.0
Δt=[t0]
Y=[y0]

for i in 1:n    
    k1 = ((((Δt[i])^2)-1)*Y[i])h
    k2 = ((((Δt[i]+(h/2))^2)-1)*(Y[i]+k1/2))h
    k3 = ((((Δt[i]+(h/2))^2)-1)*(Y[i]+k2/2))h
    k4 = ((((Δt[i]+h)^2)-1)*(Y[i]+k3))h
        
    global y0 = Y[i] + (k1+2k2+2k3+k4)/6 # updating y
    global t0 = Δt[i] + h # updating time value    
    append!(Y,y0); append!(Δt,t0) # populating arrays 
end

λ = []
for β in 1:n+1
    append!(λ,exp(-Δt[β]+(((Δt[β])^3)/3)))        
end

for j in 1:n+1
    # t, rk4-aproximation, exact-solution, absolute error
    println(Δt[j], Y[j], λ[j], abs(Y[j]-λ[j]))
end


#=
Output:
(t)  (rk4-aproximation)  (exact-solution)  (absolute error)  
0.0  1.0  1.0  0.0
0.1  0.9051391512809895  0.9051390807829696  7.049801997727201e-8
0.2  0.8209170690999957  0.8209169487181869  1.2038180885465266e-7
0.30000000000000004  0.7475158270453778  0.7475156780180909  1.4902728684429434e-7
0.4  0.6847739920729041  0.6847738325717047  1.5950119947039099e-7
0.5  0.63233681927461  0.6323366621862497  1.5708836031969042e-7
0.6  0.5897835052705566  0.5897833576128504  1.4765770617319873e-7
0.7  0.5567347179866677  0.556734581678095  1.3630857276236696e-7
0.7999999999999999  0.5329471071343952  0.5329469805876073  1.2654678793744267e-7
0.8999999999999999  0.518404336705643  0.5184042166537559  1.2005188709807157e-7
0.9999999999999999  0.5134172359588078  0.513417119032592  1.1692621582071183e-7
=#