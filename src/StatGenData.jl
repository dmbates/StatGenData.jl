using DataFrames

module StatGenData
    using DataFrames
    import Base: getindex, size

    export
      GenData,

      bedfreq

    include("GenData.jl")
end
