import BattMo

export setParameter!, setParameters!

function setParameters!(file::BattMo.JSONFile, values, experiment::String)

    if experiment=="Example"
        file.object["NegativeElectrode"]["Coating"]["ActiveMaterial"]["massFraction"] = values[1]
        file.object["PositiveElectrode"]["Coating"]["ActiveMaterial"]["massFraction"] = values[2]
    else
        println("Error: Invalid experiment!")
    end

end

function setParameters!(file::BattMo.MatlabFile, values, experiment::String)
    error("Not implemented")
end
