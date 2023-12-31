## quiets concerns of R CMD check re:no binding global varaible or function
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "GENEID",
        "PC1", "PC2", "TXNAME",
        "annotationTxId", "chr", "chr.rc",
        "compatible", "confidenceType", "endMatch",
        "eqClass", "eqClassReadCount", "equal",
        "estimates", "exon_rank", "feature",
        "geneCount", "geneId", "gene_id", "gene_sid",
        "groupVar", "group_name", "intronEnds",
        "intronStarts", "newGeneClass", "txClassDescription",
        "nobs_stored", "queryElementsOutsideMaxDist",
        "queryHits.x", "queryHits.y", "readClassId",
        "readCount", "rowMaxs", "rowMins",
        "runname", "seqlengths", "startMatch", "strand.rc",
        "subjectCount", "subjectElementsOutsideMaxDist", "sum_nobs", "txId", 
        "tx_id", "tx_name", "uniqueEndLengthQuery",
        "uniqueLengthQuery", "uniqueLengthSubject",
        "score","frame","dist","nobs","index",
        "uniqueStartLengthQuery", "value", "valueGene", "valueGeneCPM",
        "variable", "junctionMatchList", "readClass.file", 
        "start.ptm", "verbose", ".insertGaps", ".","subjectList",
        "CPM", "fullLengthCounts", "uniqueCounts",
        "annotatedEnd", "annotatedJunction", "annotatedStart", "aval", 
        "dObs", "endScore", "firstExonWidth",  "junctionEndName",
        "junctionStartName", "multi_align", "par_status",
         "spliceMotif", "startScore", "theta", "totalWidth",
        "counts","intersectWidth", "newGeneId", "newTxId", 
        "newGeneId.merge","nTx","geneReadCount",
        "startSD","endSD","readCount.posStrand","alignmentStrand","newExonId",
        "newExonId.merge","readStart","readEnd","novelGene",
        "numReads","numExons","numRCs","geneReadProp","numAstart", "numAend",
        "numTstart","numTend","tx_strand_bias", "maxTxScore", 
        "NSampleReadCount", "NSampleReadProp", "NSampleTxScore", 
        "NSampleReadCount.combined", "NSampleReadCount.new", 
        "NSampleReadProp.combined", "NSampleReadProp.new", 
        "NSampleTxScore.combined", "NSampleTxScore.new", "maxTxScore.combined",
        "maxTxScore.new" ,"txScore" ,"row_id", "sample_id", 
        "sample_name", "readCount_tmp", "group", "rowID", 
        "sumReadCount", "maxTxScore.noFit", "geneid", 
        "readIds", "readId", "anyCompatible","seqname",
        "eqClassId","K","eqClassById","rcWidth", "minRC",
        "tr_dimension","txid","cumN","gene_grp_id","txlen",
        "maxTxScore.noFit.combined","maxTxScore.noFit.new",
        "txScore.noFit","n.obs","nObs_list","K_list","txids_list",
        "anyEqual","txScore.noFit","txidTemp",
        "subjectHits.y", "txNumberFiltered"
    ))
}
