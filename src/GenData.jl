## read data from a combination of .bed, .bim and .fam files

type GenData            # Genomic data from plink-style .bed files
    snpinfo::DataFrame
    faminfo::DataFrame
    gendat::Matrix{Uint8}
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
function getindex(g::GenData, i::Integer, j::Integer)
    if !(1 <= i <= size(g.faminfo, 1) && 1 <= j <= size(g.snpinfo, 1))
        error(BoundsError)
    end
    int(0x03 & (g.gendat[div(i-1,4) + 1,j]>>(rem(i,4)<<1)))
end

function bedfreq(g::GenData)
    m,n = size(g)
    counts = zeros(Int, 4, n)
    for j in 1:n, i in 1:m counts[1 + g[i,j], j] += 1 end
    counts
end
