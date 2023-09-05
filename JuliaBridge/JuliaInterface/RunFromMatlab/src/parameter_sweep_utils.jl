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

function setParameter!(file::BattMo.MatlabFile, field::String, value::Real)
    obj = file.object
    if field == "nam-volumeFraction"
        obj["model"]["NegativeElectrode"]["ActiveMaterial"]["volumeFraction"] .= value
        obj["model"]["NegativeElectrode"]["ActiveMaterial"]["porosity"] = 1 - file.obj["model"]["PositiveElectrode"]["ActiveMaterial"]["volumeFraction"]   
    elseif field == "pam-volumeFraction"
        obj["model"]["PositiveElectrode"]["ActiveMaterial"]["volumeFraction"] .= value
        obj["model"]["PositiveElectrode"]["ActiveMaterial"]["porosity"] = 1 - file.obj["model"]["PositiveElectrode"]["ActiveMaterial"]["volumeFraction"]
    elseif field == "pam-thickness"
        #???
    elseif field == "nam-thickness"
        #???
    else
        println("ERROR: invalid flag")
    end
end

function setParameter!(file::BattMo.JSONFile, field::String, value::Real)
    obj = file.object
    if field == "nam-volumeFraction"
        obj["NegativeElectrode"]["ActiveMaterial"]["Interface"]["volumeFraction"] = value
    elseif field == "pam-volumeFraction"
        obj["PositiveElectrode"]["ActiveMaterial"]["Interface"]["volumeFraction"] = value
    elseif field == "pam-thickness"
        obj["NegativeElectrode"]["ActiveMaterial"]["thickness"] = value
    elseif field == "nam-thickness"
        obj["PositiveElectrode"]["ActiveMaterial"]["thickness"] = value
    else
        println("ERROR: invalid flag")
    end
end
