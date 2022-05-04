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

#Defino funciones
mutable struct DemoObserver <: AbstractObserver
   energy_tol::Float64
   last_energy::Float64

   DemoObserver(energy_tol=0.0) = new(energy_tol,1000.0)
end

function ITensors.checkdone!(o::DemoObserver;kwargs...)
  sw = kwargs[:sweep]
  energy = kwargs[:energy]
  if abs(energy-o.last_energy)/abs(energy) < o.energy_tol
    println("Stopping DMRG after sweep $sw")
    return true
  elseif isfile("stop") == true
    println("\n PARANDO EL DMRG POR GENERACION DEL ARCHIVO: STOP AL SWEEP $sw \n")
    return true
  end
  # Otherwise, update last_energy and keep going
  o.last_energy = energy
  return false
end

function ITensors.checkdone!(o::DMRGObserver;kwargs...)
  sw = kwargs[:sweep]

  if isfile("stop") == true
    println("\n PARANDO EL DMRG POR GENERACION DEL ARCHIVO: STOP AL SWEEP $sw \n")
    return true
  end
  return false
end


function ITensors.measure!(o::DemoObserver; kwargs...)
  energy = kwargs[:energy]
  sweep = kwargs[:sweep]
  #bond = kwargs[:bond]
  outputlevel = kwargs[:outputlevel]

  if outputlevel > 0
    #println("Sweep $sweep at bond $bond, the energy is $energy")
    println("Sweep $sweep, the energy is $energy")
  end
end

#mide valores medios
function meanvalues(N::Int64, psi::MPS, sites::Array, MVSz::Array{Float64},
                    MVNup::Array{Float64},MVNdn::Array{Float64},
                    MVNt::Array{Float64}, MVSy::Array{Float64}, MVSx::Array{Float64})

  for j=1:length(psi)
    orthogonalize!(psi, j)
    s = siteind(psi,j)
    psidag_j = dag(prime(psi[j], s))
    MVSz[j] = scalar(psidag_j * op(sites,"Sz",j) * psi[j])
    MVNup[j] = scalar(psidag_j * op(sites, "Nup", j) * psi[j])
    MVNdn[j] = scalar(psidag_j * op(sites, "Ndn", j) * psi[j])
    MVNt[j] = MVNup[j] + MVNdn[j]
    valp = scalar(psidag_j * op(sites,"Sp",j) * psi[j])
    valm = scalar(psidag_j * op(sites,"Sm",j) * psi[j])
    MVSy[j] = real((valp-valm)/(2.0im))
    MVSx[j] = real((valp+valm)/2.0)

  end
  return MVSz, MVNup, MVNdn, MVNt, MVSy, MVSx
end

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
  N = 4
  λ = 0.1
  k = 0.5


  #algoritmo para SVD ("divide_and_conquer" ,"qr_iteration", "recursive")
  # divide_and_conquer (pred) a veces tiene problema con matrices mal condicionadas
  # por eso dejo como usado a este, que es más lento
  # ver de cambiar este alg para hamiltoniano no hermitico
  alg = "qr_iteration"

  #parametros del dmrg
  sweeps = Sweeps(20)
  minsweeps = 5
  maxdim!(sweeps, 50,100,200)
  #cutoff!(sweeps, 1E-12)
  etol=1E-4  #tolerancia de energia

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


energy, psi = dmrg(H, psi0, sweeps;svd_alg=alg, outputlevel=1, ishermitian=false)
println("Final energy = $energy")


end
