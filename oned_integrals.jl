using FFTW
using Plots
using QuadGK
using MultiQuad

k_max = 400.0
n_points = 2^16+2 # 2^b + 2 so that the fft length is power of 2^(b+1)
k = LinRange(0,k_max,n_points)
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
n_fft = length(fc_t)
intensity = (0:n_fft) .* π/k_max
p = fft(fc_t)
p5 = plot(intensity[2:2:end],real.(p[2:2:end]),xaxis=("I",(0,3)),yaxis=("p(I)",(0,1000)))