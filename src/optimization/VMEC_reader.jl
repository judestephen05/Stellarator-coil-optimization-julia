using Base.RegexMatch

function read_vmec_input(input.LandremanPaul2021_QA)
    RBC = Dict{Tuple{Int, Int}, Float64}()
    ZBS = Dict{Tuple{Int, Int}, Float64}()
    NFP = 0

    for line in eachline(input.LandremanPaul2021_QA)
        if occursin("NFP =", line)
            NFP = parse(Int, split(line, "=")[2])
        end
        if occursin("MPOL =", line)
            MPOL = parse(Int, split(line, "=")[2])
        end
        if occursin("NTOR =", line)
            NTOR = parse(Int, split(line, "=")[2])
        end

        r = match(r"RBC\(\s*(-?\d+),\s*(-?\d+)\)\s*=\s*([+-]?\d+\.\d+[eE][+-]?\d+)", line)
        if r !== nothing 
            n, m, val = parse(Int, r[1]), parse(Int, r[2]), parse(Float64, r[3])
            RBC[(n, m)] = val
        end

        z = match(r"ZBS\(\s*(-?\d+),\s*(-?\d+)\)\s*=\s*([+-]?\d+\.\d+[eE][+-]?\d+)", line)
        if m !== nothing 
            n, mval, val = parse(Int, m[1]), parse(Int, m[2]), parse(Float64, m[3])
            RBC[(n, mval)] = val
        end

    end

    return NFP, RBC, ZBS
end