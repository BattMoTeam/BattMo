import JSON

export save_output

function save_output(output, outputFileName)
    stringdata = JSON.json(output)
    # write the file with the stringdata variable information
    open(outputFileName, "w") do f
        write(f, stringdata)
    end    
end