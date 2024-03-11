using Distributions, Random
using Plots

# the plan is to simulate intensity distribution by using a Poisson distribution
# over a small area.  The first thing that we have to figure out is over which length
# a normal distribution will be non-zero given double resolution

w = 2 # is the way the Gaussian width is defined in optics
sigma = w/2
ep = 1.0
x = range(-20,stop=20,length=200)
d = Normal(0.0,sigma)
intensity = ep*sigma*sqrt(2*π)*pdf.(d,x)
println(x[2]," ",x[end-1])

# set limits and define the length of the box
xlim1 = x[2]
xlim2 = x[end-1]
L = xlim2-xlim1
println(L)

# we want a certain number of particles per length
# for a particular w which represents an effective size of sqrt(pi)/2*w
N_avg = 1
N_avg_L = N_avg*L/sqrt(π)*sqrt(2)/w
print(N_avg_L)

#empty list of intensities
int_list = zeros(0)

N_samples = 1000000
p = Poisson(N_avg_L)
N_draws = rand(p,N_samples)
for i = 1:N_samples
    positions = rand(N_draws[i]).*L .+ xlim1
    intensities = ep .* sigma .* sqrt(2*π)*pdf.(d,positions)
    append!( int_list, sum(intensities) )
end

histogram(int_list,bins=500,normalize=true,xlim=(0,3),ylim=(0,1.5))

