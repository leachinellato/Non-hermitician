using Revise
using ITensors
using DelimitedFiles
using CSV
using DataFrames
using Printf
using LaTeXStrings
using Plots
using ITensors.HDF5
using MKL
using KrylovKit
using LinearAlgebra


BLAS.set_num_threads(8)

#Hamiltonian
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
  #-----------------------------------------------------------------------
  if isfile("stop")==true
      println("\n Delete stop archive!")
      stop()
  end

  #Models parameters
  N = 20
  λ = 0.1
  k = 0.5
    
  sites = siteinds("S=1/2",N;conserve_qns=false)
  

  #dmrg parameters
  alg = "qr_iteration"
  sweeps = Sweeps(1000)
  minsweeps = 5
  etol=1E-6
    
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

#obs = DMRGObserver(["Sz"], sites; energy_tol=etol,minsweeps=minsweeps)
obs = DMRGObserver(; energy_tol=etol, minsweeps=2)

energy, psi = dmrg(H, psi0, sweeps;svd_alg=alg,observer=obs, outputlevel=1, ishermitian=false)
println("Final energy = $energy")



end
