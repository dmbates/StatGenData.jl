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
    snpinfo = DataFrame(chr = int8(snp[:,1]),      # chromosome number
                        loc = int32(snp[:,4]),     # location
                        rs  = snp[:,2],
                        major = int8(snp[:,5]),
                        minor = int8(snp[:,6]))
    faminfo = DataFrame(ind = int([1:nsubj]),
                        fID = fam[:,1],
                        ID  = int(fam[:,2]),
                        pID = int(fam[:,3]),
                        mID = int(fam[:,4]),
                        sex = int8(fam[:,5]),
                        phe = int8(fam[:,6]))
    GenData(snpinfo, faminfo, bb)
end

## FIXME: Haven't handled the pairs, if any, in the trailing byte.
## Need to check which end the active bit pairs start from.
function bedfreq(g::GenData)
    b = g.gendat
    m,n = size(b)
    addnl = rem(size(g.snpinfo,1),4) # number of pairs in the trailing byte
    if bool(addnl) m -= 1 end   # m is now the number of full bytes per column
    counts = zeros(Int, 4, n)
    for j in 1:n
        for i in 1:m
            bij = b[i,j]
            for k in 0:2:6 counts[1 + (0x03 & (bij>>k)), j] += 1 end
        end
    end
    counts
end
