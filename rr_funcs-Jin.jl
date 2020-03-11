using Libdl
rrlib = Libdl.dlopen("C:\\rr\\install\\roadrunner\\bin\\roadrunner_c_api.dll")


struct RRVector
end

struct RRDoubleMatrix
end

struct RRStringArray
end

struct RRCData
end


function disableLoggingToConsole()
    ccall(dlsym(rrlib, :disableLoggingToConsole), cdecl, Bool, ())
end

function setConfigBool(key::String, value::Int64)
    ccall(dlsym(rrlib, :setConfigBool), cdecl, Int64, (Ptr{UInt8}, Cint), key, value)
end

function createRRInstance()
  val = ccall(dlsym(rrlib, :createRRInstance), cdecl, Ptr{Nothing}, ())
  if val == C_NULL
    error("Failed to start up roadRunner")
  end
  return val
end

function addCompartment(rr, cid::String, initVolume::Float64, regen::Bool)
  status = false
  if regen == true
    status = ccall(dlsym(rrlib, :addCompartment), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Cdouble), rr, cid, initVolume)
  else
    status = ccall(dlsym(rrlib, :addCompartmentNoRegen), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Cdouble), rr, cid, initVolume)
  end
  if status == false
    error(getLastError())
  end
end

function addParameter(rr::Ptr{Nothing}, pid::String, value::Float64, forceRegen::Bool)
  status = false
  if forceRegen == true
     status = ccall(dlsym(rrlib, :addParameter), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Cdouble), rr, pid, value)
  else
    status = ccall(dlsym(rrlib, :addParameterNoRegen), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Cdouble), rr, pid, value)
  end
  if status == false
    error(getLastError())
  end
end

function removeParameter(rr::Ptr{Nothing}, pid::String, forceRegen::Bool)
  status = false
  if forceRegen == true
     status = ccall(dlsym(rrlib, :removeParameter), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}), rr, pid)
  else
    status = ccall(dlsym(rrlib, :removeParameterNoRegen), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}), rr, pid)
  end
  if status == false
    error(getLastError())
  end
end

function addSpecies(rr::Ptr{Nothing}, sid::String, compartment::String, initialAmount::Float64, substanceUnit::String, regen::Bool)
  status = false
  if regen == true
    status = ccall(Libdl.dlsym(rrlib, :addSpecies), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Ptr{UInt8}, Cdouble, Ptr{UInt8}), rr, sid, compartment, initialAmount, substanceUnit)
  else
    status = ccall(Libdl.dlsym(rrlib, :addSpeciesNoRegen), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Ptr{UInt8}, Cdouble, Ptr{UInt8}), rr, sid, compartment, initialAmount, substanceUnit)
  end
  if status == false
    error(getLastError())
  end
end

function addReaction(rr::Ptr{Nothing}, rid::String, reactants::Array{String}, products::Array{String}, kineticLaw::String, regen::Bool)
  numReactants = length(reactants)
  numProducts = length(products)
  status = false
  if regen == true
    status = ccall(dlsym(rrlib, :addReaction), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Ptr{Ptr{UInt8}}, Cint, Ptr{Ptr{UInt8}}, Cint, Ptr{UInt8}), rr, rid, reactants, numReactants, products, numProducts, kineticLaw)
  else
    status = ccall(dlsym(rrlib, :addReactionNoRegen), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Ptr{Ptr{UInt8}}, Cint, Ptr{Ptr{UInt8}}, Cint, Ptr{UInt8}), rr, rid, reactants, numReactants, products, numProducts, kineticLaw)
  end
  if status == false
    error(getLastError())
  end
end

function removeReaction(rr::Ptr{Nothing}, rid::String, regen::Bool)
  status = false
  if regen == true
    status = ccall(dlsym(rrlib, :removeReaction), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}), rr, rid)
  else
    status = ccall(dlsym(rrlib, :removeReactionNoRegen), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}), rr, rid)
  end
  if status == false
    error(getLastError())
  end
end

function getNumberOfReactions(rr::Ptr{Nothing})
  return ccall(dlsym(rrlib, :getNumberOfReactions), cdecl, Int64, (Ptr{Nothing},), rr)
end

function getReactionIds(rr)
  return ccall(dlsym(rrlib, :getReactionIds), cdecl, Ptr{RRStringArray}, (Ptr{Nothing},), rr)
end

function getSBML(rr)
  char_pointer=ccall(dlsym(rrlib, :getSBML), cdecl, Ptr{UInt8}, (Ptr{Nothing},), rr)
  julia_str=unsafe_string(char_pointer)
  freeText(char_pointer)
  return julia_str
end

# function getSBML(rr)
#   return unsafe_string(ccall(dlsym(rrlib, :getSBML), cdecl, Ptr{UInt8}, (Ptr{Nothing},), rr))
# end


function getCurrentSBML(rr)
  char_pointer=ccall(dlsym(rrlib, :getCurrentSBML), cdecl, Ptr{UInt8}, (Ptr{Nothing},), rr)
  julia_str=unsafe_string(char_pointer)
  freeText(char_pointer)
  return julia_str
end

# function getCurrentSBML(rr)
#   return unsafe_string(ccall(dlsym(rrlib, :getCurrentSBML), cdecl, Ptr{UInt8}, (Ptr{Nothing},), rr))
# end

function loadSBML(rr::Ptr{Nothing}, sbml::String)
  return ccall(dlsym(rrlib, :loadSBML), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}), rr, sbml)
end

function getNumberOfGlobalParameters(rr)
  return ccall(dlsym(rrlib, :getNumberOfGlobalParameters), cdecl, Int64, (Ptr{Nothing},), rr)
end

function getNumberOfFloatingSpecies(rr)
  return ccall(dlsym(rrlib, :getNumberOfFloatingSpecies), cdecl, Int64, (Ptr{Nothing},), rr)
end

function setFloatingSpeciesInitialConcentrationByIndex(rr, index::Int64, value::Float64)
  return ccall(dlsym(rrlib, :setFloatingSpeciesInitialConcentrationByIndex), cdecl, Bool, (Ptr{Nothing}, Int64, Float64), rr, index, value)
end

function setFloatingSpeciesByIndex(rr::Ptr{Nothing}, index::Int64, value::Float64)
  status = ccall(dlsym(rrlib, :setFloatingSpeciesByIndex), cdecl, Bool, (Ptr{Nothing}, Int64, Float64), rr, index, value)
  if status == false
    error(getLastError())
  end
end

function setGlobalParameterByIndex(rr, index::Int64, value::Float64)
  return ccall(dlsym(rrlib, :setGlobalParameterByIndex), cdecl, Bool, (Ptr{Nothing}, Int64, Float64), rr, index, value)
end

function setConfigInt(key::String, value::Int64)
    ccall(dlsym(rrlib, :setConfigInt), cdecl, Int64, (Ptr{UInt8}, Cint), key, value)
end

function steadyState(rr)
  value = Array{Float64}(undef, 1)
  status = ccall(dlsym(rrlib, :steadyState), cdecl, Bool, (Ptr{Nothing}, Ptr{Float64}), rr, value)
  if status == false
    error(getLastError())
  end
  return value[1]
end

function getFloatingSpeciesConcentrations(rr)
  return ccall(dlsym(rrlib, :getFloatingSpeciesConcentrations), cdecl, Ptr{RRVector}, (Ptr{Nothing},), rr)
end

function getVectorElement(vector::Ptr{RRVector}, index::Int64)
  value = Array{Float64}(undef,1)
  ccall(dlsym(rrlib, :getVectorElement), cdecl, Bool, (Ptr{RRVector}, Int64, Ptr{Cdouble}), vector, index, value)
  return value[1]
end


function freeVector(vector::Ptr{RRVector})
  status = ccall(dlsym(rrlib, :freeVector), cdecl, Bool, (Ptr{RRVector},), vector)
  if status == false
    (error(getLastError()))
  end
end

function getScaledConcentrationControlCoefficientMatrix(rr)
  return ccall(dlsym(rrlib, :getScaledConcentrationControlCoefficientMatrix), cdecl, Ptr{RRDoubleMatrix}, (Ptr{Nothing},), rr)
end

function getMatrixElement(m::Ptr{RRDoubleMatrix}, r::Int64, c::Int64)
  value = Array{Float64}(undef,1)
  status = ccall(dlsym(rrlib, :getMatrixElement), cdecl, Bool, (Ptr{RRDoubleMatrix}, Int64, Int64, Ptr{Cdouble}), m, r, c, value)
  if status == false
    error(getLastError())
  end
  return value[1]
end

function setMatrixElement(m::Ptr{RRDoubleMatrix}, r::Int64, c::Int64, value::Float64)
  status = ccall(dlsym(rrlib, :setMatrixElement), cdecl, Bool, (Ptr{RRDoubleMatrix}, Int64, Int64, Float64), m, r, c, value)
  if status == false
    error(getLastError())
  end
end

function getGlobalParameterIds(rr)
  return ccall(dlsym(rrlib, :getGlobalParameterIds), cdecl, Ptr{RRStringArray}, (Ptr{Nothing},), rr)
end


function getStringElement(list::Ptr{RRStringArray}, index::Int64)
  char_pointer=ccall(dlsym(rrlib, :getStringElement), cdecl, Ptr{UInt8}, (Ptr{RRStringArray}, Int64), list, index)
  julia_str=unsafe_string(char_pointer)
  freeText(char_pointer)
  return julia_str
end

# function getStringElement(list::Ptr{RRStringArray}, index::Int64)
#   return unsafe_string(ccall(dlsym(rrlib, :getStringElement), cdecl, Ptr{UInt8}, (Ptr{RRStringArray}, Int64), list, index))
# end


function freeRRInstance(rr::Ptr{Nothing})
  free_status = ccall(dlsym(rrlib, :freeRRInstance), cdecl, Bool, (Ptr{Nothing},), rr)
  if free_status == false
    error(getLastError())
  end
end

function freeMatrix(matrix::Ptr{RRDoubleMatrix})
  status = ccall(dlsym(rrlib, :freeMatrix), cdecl, Bool, (Ptr{RRDoubleMatrix},), matrix)
  if status == false
    (error(getLastError()))
  end
end

# does not work, why?
# function getLastError()
#   char_pointer=ccall(dlsym(rrlib, :getLastError), cdecl, Ptr{UInt8}, ())
#   julia_str=unsafe_string(char_pointer)
#   freeText(char_pointer)
#   return julia_str
# end


function getLastError()
  return unsafe_string(ccall(dlsym(rrlib, :getLastError), cdecl, Ptr{UInt8}, ()))
end

function simulate(rr::Ptr{Nothing})
  return ccall(dlsym(rrlib, :simulate), cdecl, Ptr{RRCData}, (Ptr{Nothing},), rr)
end

function getRateOfChange(rr, index::Int64)
  value = Array{Float64}(undef,1)
  ccall(dlsym(rrlib, :getRateOfChange), cdecl, Bool, (Ptr{Nothing}, Int64, Ptr{Float64}), rr, index, value)
  return value[1]
end

function freeRRCData(handle::Ptr{RRCData})
  status = ccall(dlsym(rrlib, :freeRRCData), cdecl, Bool, (Ptr{RRCData},), handle)
  if status == false
    error(getLastError())
  end
end

function freeText(text::Ptr{UInt8})
  status = ccall(dlsym(rrlib, :freeText), cdecl, Bool, (Ptr{UInt8},), text)
  if status == false
    error(getLastError())
  end
end

function freeStringArray(sl::Ptr{RRStringArray})
  status =  ccall(dlsym(rrlib, :freeStringArray), cdecl, Bool, (Ptr{RRStringArray},), sl)
  if status == false
    error(getLastError())
  end
end

function getFloatingSpeciesIds(rr)
  return ccall(dlsym(rrlib, :getFloatingSpeciesIds), cdecl, Ptr{RRStringArray}, (Ptr{Nothing},), rr)
end


function getBoundarySpeciesIds(rr)
  return ccall(dlsym(rrlib, :getBoundarySpeciesIds), cdecl, Ptr{RRStringArray}, (Ptr{Nothing},), rr)
end

function getNumberOfBoundarySpecies(rr)
  return ccall(dlsym(rrlib, :getNumberOfBoundarySpecies), cdecl, Int64, (Ptr{Nothing},), rr)
end

function getBoundarySpeciesConcentrations(rr)
  return ccall(dlsym(rrlib, :getBoundarySpeciesConcentrations), cdecl, Ptr{RRVector}, (Ptr{Nothing},), rr)
end

function setBoundarySpeciesByIndex(rr::Ptr{Nothing}, index::Int64, value::Float64)
  status = ccall(dlsym(rrlib, :setBoundarySpeciesByIndex), cdecl, Bool, (Ptr{Nothing}, Int64, Float64), rr, index, value)
  if status == false
    error(getLastError())
  end
end

function resetRR(rr)
  return ccall(dlsym(rrlib, :reset), cdecl, Bool, (Ptr{Nothing},), rr)
end

function resetAllRR(rr::Ptr{Nothing})
  return ccall(dlsym(rrlib, :resetAll), cdecl, Bool, (Ptr{Nothing},), rr)
  if status == false
    error(getLastError())
  end
end

function setBoundary(rr::Ptr{Nothing}, sid::String, boundaryCondition::Bool, forceRegen::Bool)
  status = false
  if forceRegen == true
    status = ccall(dlsym(rrlib, :setBoundary), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Bool), rr, sid, boundaryCondition)
  else
    status = ccall(dlsym(rrlib, :setBoundaryNoRegen), cdecl, Bool, (Ptr{Nothing}, Ptr{UInt8}, Bool), rr, sid, boundaryCondition)
  end
  if status == false
    error(getLastError())
  end
end
