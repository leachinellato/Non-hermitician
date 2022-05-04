using Revise
using ITensors
using Printf
using MKL
using KrylovKit
using LinearAlgebra


BLAS.set_num_threads(8)

#define the hamiltonian

function Hamiltonian(sites, λ::Float64, k::Float64)

    N = length(sites)
    ampo = AutoMPO()

    for j in 1:(N - 1)
        ampo += -0.5*2*λ,"Sz",j,"Sz",j+1
    end
    for j in 1:N
        ampo +=  -0.5*2*im*k, "Sz",j
        ampo += -0.5, "S+",j
        ampo += -0.5, "S-",j
    end
  return MPO(ampo, sites)
end




let

  #models parameters
  N = 4
  λ = 0.1
  k = 0.5

  alg = "qr_iteration"

  # DMRG parameters
  sweeps = Sweeps(20)
  minsweeps = 5
  maxdim!(sweeps, 50,100,200)
  #cutoff!(sweeps, 1E-12)
  etol=1E-4 

#-----------------------------------------------------------------------
 sites = siteinds("S=1/2",N;conserve_qns=false)
#-----------------------------------------------------------------------
   
#Estado inicial
state = ["Emp" for n in 1:N]
p = N
for i in N:-1:1
  if p > i
    #println("Doubly occupying site $i")
    state[i] = "UpDn"
    p -= 2
elseif p > 0
    #println("Singly occupying site $i")
    state[i] = (isodd(i) ? "Up" : "Dn")
    p -= 1
  end
end
psi0 = randomMPS(sites, state)
@show flux(psi0)
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

H = Hamiltonian(sites, λ, k)
#-----------------------------------------------------------------------

energy, psi = dmrg(H, psi0, sweeps;svd_alg=alg, outputlevel=1, ishermitian=false)
println("Final energy = $energy")


end
