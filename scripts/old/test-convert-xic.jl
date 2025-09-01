
using OrderedCollections
using DelimitedFiles
using YAML

description = readdlm(joinpath(@__DIR__, "..", "data", "raw-xic-values.dlm"))

function consistent_label(name, coupling)
    if startswith(name, "K*")
        coupling == "1/2,0" && return 1
        coupling == "1/2,-1" && return 2
        coupling == "-1/2,1" && return 3
        coupling == "-1/2,0" && return 4
    end
    if startswith(name, "K0")
        coupling == "1/2,0" && return 2
        coupling == "-1/2,0" && return 1
    end
    if startswith(name, "Λ") || startswith(name, "∆")
        coupling == "1/2,0" && return 1
        coupling == "-1/2,0" && return 2
    end
    return -1
end

processed_description = map(eachrow(description)) do row
    type, name, couplings = row[1], row[2], row[3]
    if last(type) == 'H'
        index = consistent_label(name, couplings)
        _reim = (row[1] == "ReH") ? "r" : ((row[1] == "ImH") ? "i" : "")
        _name = replace(name, "Λ" => "L", "∆" => "D", "*" => "", "+" => "", ")0" => ")", "0(" => "(")
        return "A$(_reim)$(_name)$(index)" => "$(row[4]) ± $(row[5]) ± $(row[6]) ± $(row[7])"
    end
    if first(type) == 'P'
        _name = replace(name, "γ" => "gamma", "*" => "", "+" => "", ")0" => ")", "0(" => "(")
        return "$(_name)" => "$(row[4]) ± $(row[5]) ± $(row[6]) ± $(row[7])"
    end
    error("Unknown name format: $type")
end

YAML.write_file(joinpath(@__DIR__, "..", "data", "xic2pKpi.yaml"), Dict("parameters" => LittleDict(processed_description...)))
