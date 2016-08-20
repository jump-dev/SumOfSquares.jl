# Inspired from SOSTools
import Base.extrema

function extrema(f, itr)
  state = start(itr)
  it1, state = next(itr, state)
  mini = maxi = f(it1)
  while !done(itr, state)
    it, state = next(itr, state)
    fit = f(it)
    if fit < mini
      mini = fit
    end
    if fit > maxi
      maxi = fit
    end
  end
  mini, maxi
end

cfld(x::NTuple{2,Int}, n) = (cld(x[1], 2), fld(x[2], 2))

# TODO sparse with Newton polytope (Polyhedra.jl for convex hull)
function getmonomialsforcertificate(Z::MonomialVector, sparse=:No)
  if sparse == :No
    # Cheap approximation of the convex hull as the approximation of:
    #
    # z such that mindeg < sum(z) < maxdeg
    # |\
    # |#\ <-------- sum(z) = maxdeg
    # |##\
    # |\<-\-------- sum(z) = mindeg
    # | \##\
    # +---------
    #
    # and:
    #
    # z such that minmultideg < z < maxmultideg
    # | +----+ <--- maxmultideg
    # | |####|
    # | |####|
    # | +----+
    # | ^---------- minmultideg
    # +---------
    mindeg, maxdeg = cfld(extrema(map(sum, Z.Z)), 2)
    n = length(Z.Z[1])
    minmultideg, maxmultideg = Vector{Int}(n), Vector{Int}(n)
    for i in 1:n
      minmultideg[i], maxmultideg[i] = cfld(extrema(z->z[i], Z.Z), 2)
    end
    MonomialVector(vars(Z), mindeg:maxdeg, z -> reduce(&, true, minmultideg .<= z .<= maxmultideg))
  else
    error("Not supported yet :(")
  end
end
getmonomialsforcertificate(Z::Vector, sparse=:No) = getmonomialsforcertificate(MonomialVector(Z), sparse)
