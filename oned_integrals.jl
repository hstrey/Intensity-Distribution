using FFTW
using CairoMakie
using QuadGK
using MultiQuad
using Distributions, Random

k_max = 400.0
n_points = 2^16+2 # 2^b + 2 so that the fft length is power of 2^(b+1)

function oned_cfarg(k_max,npoints)
    k = LinRange(0,k_max,n_points)
    f_real = Float64[]
    f_imag = Float64[]
    for kk in k
        value, _ = quadgk(x->(cos(kk*x)-1)/x/sqrt(-log(x)),0,1.0)
        append!(f_real,value)
        value, _ = quad(x->sin(kk*x)/x/sqrt(-log(x)),0,1.0,method=:gausslegendre, order=10000)
        append!(f_imag,value)
    end
    return f_real, f_imag
end

f_real, f_imag = oned_cfarg(400.0,2^16+2)

function oned_intprob(f_reaL, f_imag, w, rho)
    # prefactor is 1/sqrt(2)*w*c for 1-dim case
    prefactor = rho*w/sqrt(2)
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
    p = fft(fc_t)/n_fft*400/π*2 # normalized fft
    return intensity, p
end

# the plan is to simulate intensity distribution by using a Poisson distribution
# over a small area.  The first thing that we have to figure out is over which length
# a normal distribution will be non-zero given double resolution
function oned_simul(xlim1,xlim2,w,rho,samples)

    sigma = w/2
    ep = 1.0
    d = Normal(0.0,sigma)

    # set limits and define the length of the box
    L = xlim2-xlim1

    # we want a certain number of particles per length
    # for a particular w which represents an effective size of sqrt(pi)/2*w
    N_avg_L = rho*L

    #empty list of intensities
    int_list = Float64[]

    pd = Poisson(N_avg_L)
    N_draws = rand(pd,samples)
    for i = 1:samples
        positions = rand(N_draws[i]).*L .+ xlim1
        intensities = ep .* sigma .* sqrt(2*π)*pdf.(d,positions)
        append!( int_list, sum(intensities) )
    end
    return int_list
end

w = 1.0 # is the way the Gaussian width is defined in optics
rho = 0.2 # density of particles N/L
intensity, p = oned_intprob(f_real, f_imag, w, rho)

# set limits and define the length of the box
x = range(-20,stop=20,length=200)
xlim1 = x[2]
xlim2 = x[end-1]

int_list = oned_simul(xlim1,xlim2,w,rho,1000000)

function distfig(int_list,intensity,p)
    f = Figure()
    ax = Axis(f[1, 1],
        title = "1-d Intensity Distribution",
        xlabel = "I",
        ylabel = "p(I)"
        )
    xlims!(ax,0,3)
    ylims!(ax,0,2)
    hist!(ax,int_list,
            bins=200,
            strokewidth = 1, strokecolor = :darkgrey,
            color=:gray,
            normalization = :pdf,
            label="simulation")
    lines!(ax,intensity[2:2:end],
            real.(p[2:2:end]),
            color=:black,
            linewidth=2,
            linestyle=:solid,
            label="theory")
    axislegend(ax)
    text!(1.5, 1; text="ρ=0.2",fontsize = 18)
    f
end
f = distfig(int_list, intensity, p)
save("oned_rho02Makie.png",f)
