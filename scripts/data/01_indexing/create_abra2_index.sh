mkdir -p ../../../test/data/indices/abra2/
nexus_create_abra2_targets_bed_file \
  --gencode_gtf_file ../../../test/data/gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf \
  --bedtools bedtools \
  --output_temp_bed_file ../../../test/data/indices/abra2/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.abra2_targets_temp.bed \
  --output_bed_file ../../../test/data/indices/abra2/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.abra2_targets.bed