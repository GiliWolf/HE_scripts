
file_seperator="_"
mate_seperator="-"
suffix_seperator="."
PE=1
mate1="A2C_ADAR_GMCSF_AdarWT_MDA5KO_79_372_SKO-2.fq"

if [ "$PE" -eq 0 ]; then
    sample_id=$(echo "${mate1}" | awk -v fs="${file_seperator}" -v ss="${suffix_seperator}" 'BEGIN{FS=ss}{split($1, parts, fs); for(i=2; i<length(parts); i++) printf "%s_", parts[i]; printf "%s", parts[length(parts)]}')
else
    if [[ "$mate_seperator" == "$file_seperator" ]]; then
        sample_id=$(echo "${mate1}" | awk -v fs="${file_seperator}" -v ss="${suffix_seperator}" 'BEGIN{FS=ss}{split($1, parts, fs); for(i=2; i<(length(parts)-1); i++) printf "%s_", parts[i]; printf "%s", parts[length(parts)-1]}')
    else
        sample_id=$(echo "${mate1}" | awk -v fs="${file_seperator}" -v ms="${mate_seperator}" -v ss="${suffix_seperator}" 'BEGIN{FS=ms ss}{split($1, parts, fs); for(i=2; i<length(parts); i++) printf "%s_", parts[i]; printf "%s", parts[length(parts)]}')
    fi

fi

echo "$sample_id"




