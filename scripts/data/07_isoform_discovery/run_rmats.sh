nexus run --nf-workflow isoform_characterization_rmats.nf \
  -c /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/nextflow/nextflow_test_docker.config \
  -w /Users/leework/Documents/Research/projects/project_nexus/data/processed/work/isoform_discovery_rmats/ \
  --samples_tsv_file /Users/leework/Documents/Research/projects/project_nexus/nexus/scripts/data/07_isoform_discovery/samples/samples_rmats.tsv \
  --reference_genes_gtf_file /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/gtf/gencode.v45.chr_patch_hapl_scaff.annotation.chr17.gtf.gz \
  --params_rmats '"-t paired --libType fr-unstranded --readLength 151"' \
  --output_dir /Users/leework/Documents/Research/projects/project_nexus/nexus/test/data/inputs/mopepgen/
