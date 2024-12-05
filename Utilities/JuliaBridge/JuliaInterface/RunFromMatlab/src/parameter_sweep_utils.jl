import BattMo

export setParameter!, setParameters!

function setParameters!(input::BattMo.InputParams, values, experiment::String)

    if experiment=="Example"
        input["NegativeElectrode"]["Coating"]["ActiveMaterial"]["massFraction"] = values[1]
        input["PositiveElectrode"]["Coating"]["ActiveMaterial"]["massFraction"] = values[2]
    else
        println("Error: Invalid experiment!")
    end

end

function setParameters!(file::BattMo.MatlabInputParams, values, experiment::String)
    error("Not implemented")
end
