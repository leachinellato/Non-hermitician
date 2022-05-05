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



#define el hamiltoniano

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
#  Convert these terms to an MPO tensor network
  return MPO(ampo, sites)
end




let
  #-----------------------------------------------------------------------
  #Chequea que no exista el archivo stop, el cual mata la corrida del DMRG
  if isfile("stop")==true
      println("\n Borrar el archivo stop!")
      stop()
  end

  #parametros del modelo
  N = 20
  λ = 0.1
  k = 0.5


  #algoritmo para SVD ("divide_and_conquer" ,"qr_iteration", "recursive")
  # divide_and_conquer (pred) a veces tiene problema con matrices mal condicionadas
  # por eso dejo como usado a este, que es más lento
  # ver de cambiar este alg para hamiltoniano no hermitico
  alg = "qr_iteration"

  #parametros del dmrg
  sweeps = Sweeps(1000)
  minsweeps = 5
  maxdim!(sweeps, 50,100,200)
  #cutoff!(sweeps, 1E-12)
  etol=1E-12  #tolerancia de energia

  #-----------------------------------------------------------------------


 sites = siteinds("S=1/2",N;conserve_qns=false)
#-----------------------------------------------------------------------
#Observadores de energía
#obs = DMRGObserver(["Sz"], sites; energy_tol=etol,minsweeps=minsweeps)
#obs = DemoObserver(etol)  #defino el operador observador de energia
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

#obs = DMRGObserver(["Sz"], sites; energy_tol=etol,minsweeps=minsweeps)
obs = DMRGObserver(; energy_tol=etol, minsweeps=2)

energy, psi = dmrg(H, psi0, sweeps;svd_alg=alg,observer=obs, outputlevel=1, ishermitian=false)
println("Final energy = $energy")



end
