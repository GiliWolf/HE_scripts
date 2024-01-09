# output files:
ref_base='A'
alt_base='G'
fastq_path="/home/alu/fulther/Validations/GiliTest/first_part/pre_analysis/unmapped/SRR11548778/Unmapped.out.mate1"
fastq_id="SRR11548778_1_Unmapped"

transformed_fastq_path="./${fastq_id}_${ref_base}2${alt_base}.fastq"

        # transform each letter in the sequences of a fastq file from ref_bas > alt_base
        transform_fastq() {
            #create lowercase ref and alt bases
            ref_base_lower=$(echo "${ref_base}" | tr '[:upper:]' '[:lower:]')
            alt_base_lower=$(echo "${alt_base}" | tr '[:upper:]' '[:lower:]')
            #create uppercase ref and alt bases
            ref_base_upper=$(echo "${ref_base}" | tr '[:lower:]' '[:upper:]')
            alt_base_upper=$(echo "${alt_base}" | tr '[:lower:]' '[:upper:]')

            # modify only sequence records (line) using the format 4k+2 (start on line 2, repeat each 4 lines)
            # transform all lowercase ref bases into lowercase alt bases
            # transform all uppercase ref bases into uppercase alt bases
            awk -v ref_base_lower="${ref_base_lower}" -v alt_base_lower="${alt_base_lower}" \
            -v ref_base_upper="${ref_base_upper}" -v alt_base_upper="${alt_base_upper}" \
            '
            NR%4==2 {
                # Transform lowercase ref bases into lowercase alt bases
                gsub(ref_base_lower, alt_base_lower);
                # Transform uppercase ref bases into uppercase alt bases
                gsub(ref_base_upper, alt_base_upper)
            } 1' "${fastq_path}" > "${transformed_fastq_path}"
        
        }

        #transformation
            echo -e "---------------------\nstarted tranforming ${fastq_path}\n"
            transform_fastq
            echo -e "********************\nended transforming ${transformed_fastq_path}\n---------------------\n"
