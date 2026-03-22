import JSON
import Jutul
import BattMo

######################################################################################################
#JSON methods: When we go from Julia to Matlab we write the output to json file, which are then read
#by matlab. We therefore have to design custom methods for how JSON should serialize custom Types
#defined in BattMo and Jutul.
######################################################################################################

JSON.lower(c::BattMo.SimpleControllerCV) = Dict(:time              => c.time,
                                                :target            => c.target,
                                                :target_is_voltage => c.target_is_voltage)

JSON.lower(config::Jutul.JutulConfig)=Dict(:name    => config.name,
                                           :values  => config.values,
                                           :options => config.options)

JSON.lower(sim::Jutul.JutulSimulator)           = typeof(sim)
JSON.lower(obj::Nothing)                        = ""
JSON.lower(ptr::Ptr)                            = typeof(ptr)
JSON.lower(relax::Jutul.NonLinearRelaxation)    = typeof(relax)
JSON.lower(mapp::Jutul.AbstractGlobalMap)       = typeof(mapp)
JSON.lower(mapp::Jutul.TrivialGlobalMap)        = typeof(mapp)
JSON.lower(formulation::Jutul.JutulFormulation) = typeof(formulation)
JSON.lower(entity::Jutul.JutulEntity)           = typeof(entity)
JSON.lower(scalar::Jutul.ScalarVariable)        = typeof(scalar)
JSON.lower(vec::Jutul.VectorVariables)          = typeof(vec)
JSON.lower(flux::Jutul.FluxType)                = typeof(flux)
JSON.lower(eq::Jutul.JutulEquation)             = typeof(eq)
JSON.lower(domain::Jutul.JutulDomain)           = typeof(domain)
JSON.lower(sys::Jutul.JutulSystem)              = typeof(sys)

JSON.lower(func::typeof(BattMo.computeOCP_Graphite_Torchio))                     = Base.nameof(func)
JSON.lower(func::typeof(BattMo.compute_reaction_rate_constant_graphite)) = Base.nameof(func)
JSON.lower(func::typeof(BattMo.computeOCP_Graphite_SiOx_Chen2020))       = Base.nameof(func)
JSON.lower(func::typeof(BattMo.computeOCP_NMC811_Chen2020))              = Base.nameof(func)
JSON.lower(func::typeof(BattMo.computeOCP_LCO))                       = Base.nameof(func)
JSON.lower(func::typeof(BattMo.compute_reaction_rate_constant))          = Base.nameof(func)
JSON.lower(func::typeof(BattMo.computeElectrolyteConductivity_default))  = Base.nameof(func)
JSON.lower(func::typeof(BattMo.computeElectrolyteConductivity_Chen2020)) = Base.nameof(func)
JSON.lower(func::typeof(BattMo.computeDiffusionCoefficient_default))     = Base.nameof(func)
JSON.lower(func::typeof(BattMo.computeDiffusionCoefficient_Chen2020))    = Base.nameof(func)

julia_version() = string(Base.VERSION)
