using FFTW
using SpecialFunctions
using Plots

# here we will calculate the probability distribution
# from the characterisitic SpecialFunctions
k = LinRange(0,200,10000)
f_real = cosint.(k)
f_imag = sinint.(k)
p1 = plot(k,f_real)
p2 = plot!(k,f_imag)

prefactor = 0.2
er = exp(-prefactor*Base.MathConstants.eulergamma).*
     exp.(prefactor*f_real).*
     k.^(-prefactor)
ei1 = cos.(prefactor .* f_imag)
ei2 = sin.(prefactor .* f_imag)
er[1]=1
# plot er
p3 = plot(k[2:end],er[2:end],xaxis=:log, yaxis=:log)
# plot real and imaginary part
plot(k[2:end],ei1[2:end],xaxis=:log, yaxis=:log)
p4 = plot!(k[2:end],ei2[2:end],xaxis=:log, yaxis=:log)

fc1 = er .* ei1 + im * (er .* ei2)
fc2 = reverse(er .* ei1 - im * (er .* ei2))
padding = zeros(45537).*(0.0 + im*0)
fc = [fc1;fc2[1:end-1];padding]
println(length(fc1)," ",length(fc))
p = fft(fc)
p5 = plot(real.(p),xaxis=("x",(0,1000)),yaxis=("p(x)",(0,1000),:log))
theta = 2*π/65536*vcat(0:65535)
w = 2.0 .* (1.0 .- cos.(theta))/theta.^2
alpha0 = -(1.0 .- cos.(theta))/theta.^2 + im .* (theta .- sin.(theta))/theta.^2
intensity = w .* p + fc1[1] .* alpha0 + exp.(2*π*im*19999 .* vcat(0:65535))
p6 = plot(real.(intensity),xaxis=("x",(0,1000)),yaxis=("p(x)",(0,1000),:log))
println("hello")