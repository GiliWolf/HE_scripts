

def getSampleID(file_name) {
    def file_seperator="_"
    def mate_seperator="-"
    def suffix_seperator = "."
    def PE = 0
    def tokens = file_name.tokenize(suffix_seperator).get(0)
    if (PE == 1){
        // PE
        if (mate_seperator != file_seperator){
            return tokens.tokenize(mate_seperator).get(0)
        }
        else {
            def without_suffix = file_name.tokenize(mate_seperator)
            return without_suffix[0..-2].join(file_seperator)
        }
        
    }
    // SE
    else{
        return tokens
    }

}

def sampleID1  = getSampleID("ADAR_GMCSF_AdarWT_MDA5KO_79_372_SKO-1.fq")

println(sampleID1)