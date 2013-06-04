## read data from a combination of .bed, .bim and .fam files

type GenData            # Genomic data from plink-style .bed files
    snpinfo::DataFrame
    faminfo::DataFrame
    gendat::Matrix{Uint8}
end

type GenData2
    gendat::Matrix{Uint8}
    nsubj::Int
    counts::Matrix{Int}
end

function colrange(c::Vector, mv::Float64=-Inf)
    mx = -Inf; mn = Inf; iv = true; str = false; n = length(c); msng = falses(n)
    for i in 1:n
        v = c[i]
        if v == mv msng[i] = true; continue end
        if typeof(v) != Float64 str = true; continue end
        if !isinteger(v) iv = false end
        if v > mx mx = v end
        if v < mn mn = v end
    end
    str, iv, mn, mx, msng
end

function vals(c::Vector, mv::Float64=-Inf)
    str, iv, mn, mx, msng = colrange(c,mv); n = length(c)
    if str
        cc = similar(c,ASCIIString)
        for i in 1:n
            if msng[i] cc[i] = ""; continue end
            v = c[i]
            if typeof(v) <: Real cc[i] = string(iv?int(v):v); continue end
            cc[i] = string(v)
        end
        uu = unique(cc)
        if length(uu) > div(n + 3,4) return DataArray(cc,msng) end
        return compact(PooledDataArray(cc, uu, msng, Uint64))
    end
    if iv
        if mn >= 0
            if mx < typemax(Uint8) return DataArray(uint8(c), msng) end
            if mx < typemax(Uint16) return DataArray(uint16(c), msng) end
            if mx < typemax(Uint32) return DataArray(uint32(c), msng) end
            return DataArray(int(c), msng)
        end
        if typemin(Int8) <= mn <= mx <= typemax(Int8) return DataArray(int8(c), msng) end
        if typemin(Int16) <= mn <= mx <= typemax(Int16) return DataArray(int16(c), msng) end
        if typemin(Int32) <= mn <= mx <= typemax(Int32) return DataArray(int32(c), msng) end
        return DataArray(int(c), msng)
    end
    DataArray(float64(c), msng)
end

function GenData(bnm::ASCIIString) # bnm = base file name without extension
    bimnm = string(bnm,".bim")
    snp = readdlm(bimnm)
    if size(snp,2) != 6
        error("file $binnm should have 6 tab-delimited columns")
    end
    nsnp = size(snp,1)
    famnm = string(bnm,".fam")
    fam = readdlm(famnm)
    if size(fam,2) != 6
        error("file $famnm should have 6 tab-delimited columns")
    end
    nsubj = size(fam,1)
    bednm = string(bnm,".bed")
    s = open(bednm)
    b1 = read(s,Uint8); b2 = read(s,Uint8); b3 = read(s,Uint8)
    if b1 != 0x6c || b2 != 0x1b error("wrong magic number in file $bednm") end
    if b3 != 1 error(".bed file, $bednm, is not in correct orientation") end
    bb = mmap_array(Uint8, (div(nsubj+3,4),nsnp), s)
    # snpinfo is a data frame with columns
    #c1=chrom,c2=loc,c3=rs,c4=MajorAllele,c5=minor allele,c6=MAdirection(1,2),c7=MAF
    #c8=Hom11count,c9=het12count,c10=hom22count,c11=#missing
    major = typeof(snp[1,5]) <: Number ? int8(snp[:,5]) :
            typeof(snp[1,5]) <: String ?
            compact(PooledDataArray(convert(Vector{typeof(snp[1,5])},snp[:,5]),
                                    trues(nsnp),Uint8)) :
            error("type of snpinfo column 5, $(typeof(snp[1,5])), is neither Numeric nor String") 
    minor = typeof(snp[1,6]) <: Number ? int8(snp[:,6]) :
            typeof(snp[1,6]) <: String ?
            compact(PooledDataArray(convert(Vector{typeof(snp[1,6])},snp[:,6]),
                                    trues(nsnp),Uint8)) :
            error("type of snpinfo column 6, $(typeof(snp[1,6])), is neither Numeric nor String")   
    loc = max(snp[:,4]) < typemax(Int) ? int(snp[:,4]) : snp[:,4]
    snpinfo = DataFrame(chr = int8(snp[:,1]),      # chromosome number
                        loc = loc,                 # location
                        rs  = snp[:,2],
                        major = major,
                        minor = minor)
    faminfo = DataFrame(ind = int([1:nsubj]),
                        fID = fam[:,1],
                        ID  = int(fam[:,2]),
                        pID = int(fam[:,3]),
                        mID = int(fam[:,4]),
                        sex = int8(fam[:,5]), # convert to PooledDataArray with levels "M" and "F"
                        phe = int8(fam[:,6])) # check for and encode missing data values
    GenData(snpinfo, faminfo, bb)
end

size(g::GenData) = size(g.faminfo,1), size(g.snpinfo,1)
getindex(g::GenData, i::Integer, j::Integer) = 0x03 & (g.gendat[div(i-1,4) + 1,j]>>(rem(i,4)<<1))

function bedfreq(g::GenData)
    m = size(g.faminfo, 1);
    bb = g.gendat; n = size(bb,2)
    counts = zeros(Int, 4, n)
    for j in 1:n, i in 1:m
        counts[1 + (0x03 & (bb[div(i-1,4)+1,j]>>(rem(i,4)<<1))), j] += 1
    end
    counts
end

function bedfreq(b::Matrix{Uint8},ns::Integer)
    m,n = size(b)
    if div(ns+3,4) != m error("ns = $ns should be in [$(4m-3),$(4m)] for size(b) = $(size(b))") end
    counts = zeros(Int, 4, n)
    bpt = convert(Ptr{Uint8},b)
    for j in 1:n
        ptj = bpt + (j-1)*m
        for i in 1:ns
            ii = 1 + (0x03 & (unsafe_load(ptj, div(i-1,4)+1)>>(rem(i,4)<<1)))
            counts[ii,j] += 1
        end
    end
    counts
end

function GenData2(bnm::ASCIIString)     # bnm = basename of .bed, .bim and .fam files
    # counting the lines in the .bim and .fam files to get the dimensions (will clean this up)
    nsnp = int(split(readall(`wc -l $(string(bnm,".bim"))`))[1])
    nsubj = int(split(readall(`wc -l $(string(bnm,".fam"))`))[1])
    bednm = string(bnm,".bed")
    s = open(bednm)
    b1 = read(s,Uint8); b2 = read(s,Uint8); b3 = read(s,Uint8)
    if b1 != 0x6c || b2 != 0x1b error("wrong magic number in file $bednm") end
    if b3 != 1 error(".bed file, $bednm, is not in correct orientation") end
    m = div(nsubj+3,4)
    bb = mmap_array(Uint8, (m,nsnp), s)
    GenData2(bb,nsubj,bedfreq(bb,nsubj)')
end
