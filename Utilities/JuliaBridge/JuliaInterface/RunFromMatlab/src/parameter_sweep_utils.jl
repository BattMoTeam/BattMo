import BattMo

export setParameter!, setParameters!

function setParameters!(file::BattMo.JSONFile, values, experiment::String)
    if experiment=="Example"
        file.object["NegativeElectrode"]["ActiveMaterial"]["Interface"]["volumeFraction"]=values[1]
        file.object["PositiveElectrode"]["ActiveMaterial"]["Interface"]["volumeFraction"]=values[2]
    else
        println("Error: Invalid experiment!")
    end
end

function setParameters!(file::BattMo.MatlabFile, values, experiment::String)
    error("Not implemented")
end
