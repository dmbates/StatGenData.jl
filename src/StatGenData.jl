using DataFrames

module StatGenData
    using DataFrames
    import Base: getindex, size

    export
      GenData,
      GenData2,

      bedfreq

    include("GenData.jl")
end
