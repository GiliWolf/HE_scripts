

// def getSampleID(file_name) {
    def file_seperator="_"
    // def seperated_file =  file_name.toString().split("${file_seperator}Unmapped")
    // return seperated_file[0]

    def mate_seperator="_"
    def suffix_seperator = "."
    def PE = 0
    // def tokens = file_name.tokenize(suffix_seperator).get(0)
    // if (PE == 1){
    //     // PE
    //     if (mate_seperator != file_seperator){
    //         return tokens.tokenize(mate_seperator).get(0)
    //     }
    //     else {
    //         def without_suffix = file_name.tokenize(mate_seperator)
    //         return without_suffix[0..-2].join(file_seperator)
    //     }
        
    // }
    // // SE
    // else{
    //     return tokens
    // }

// }
sampleID = "ADAR_GMCSF_AdarWT_MDA5KO_79_372_SKO_1.mate1"
def sampleID_parsed  = sampleID.split("[${mate_seperator}*${suffix_seperator}]")

println(sampleID_parsed)