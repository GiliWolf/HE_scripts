
samtools view -h /private10/Projects/Gili/HE_workdir/first_part/GTEX/BrainCerebellum_PE_multimappers/re-transform/A2G/A2G_GTEX-11DZ1_BrainCerebellum_GTEX-11DZ1-2926-SM-5A5KI_re-transformed.sam | grep -P "(NH:i:[2-9]|NH:i:[1-9][0-9]|^@)" - > myBAM_filtered.sam