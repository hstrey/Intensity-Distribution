using FFTW
using PyPlot
using QuadGK
using MultiQuad

k = LinRange(0,200,16386)
f_real = zeros(0)
f_imag = zeros(0)
for kk in k
    value, err = quadgk(x->(cos(kk*x)-1)/x/sqrt(-log(x)),0,1.0)
    append!(f_real,value)
    value, err = quad(x->sin(kk*x)/x/sqrt(-log(x)),0,1.0,method=:gausslegendre, order=10000)
    append!(f_imag,value)
end

prefactor = 0.5
er = exp.(prefactor*f_real)
ei1 = cos.(prefactor .* f_imag)
ei2 = sin.(prefactor .* f_imag)

er[1]=1
fc1 = er .* ei1 + im * (er .* ei2)
fc2 = conj.(fc1)

fc = (fc1[1:end-1] .+ fc2[2:end])/2
# make sure that p(x) is normalized and that the highest frequency
# has no imaginary part
fc[1]=1
fc[end]=real(fc[end])

fc_t = [fc; reverse(conj.(fc))[2:end-1]]
println(length(fc_t))
p = fft(fc_t)
p5 = plot(real.(p),xaxis=("x",(0,200)),yaxis=("p(x)",(0,500)))