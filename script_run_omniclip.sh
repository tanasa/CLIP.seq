# 1. Generate the database

docker run --rm -it \
  -v "/home/tanasa/16_hours/samples_16h_nf_P13_peaks/COMBINED_bam_analysis":/data \
  kalasnty/omniclip_container \
  generateDB \
  --gff-file /data/gencode.v48.basic.annotation.gff3 \
  --db-file /data/gencode.v48.basic.annotation.db

# 2. Parsing the BG file

mkdir A27_13.pool_final.NM.omniclip

docker run --rm -it \
  -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/pooled_omniclip":/data \
  kalasnty/omniclip_container \
  parsingBG \
  --db-file /data/gencode.v48.basic.annotation.db \
  --bg-files /data/A27_13_In.pool_final.NM.bam \
  --out-file /data/A27_13_In.pool_final.NM.bam.out \
  --genome-dir /data/hg38.ucsc.fasta.chr


# 3. Parsing the CLIP file

docker run --rm -it \
  -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/pooled_omniclip":/data \
  kalasnty/omniclip_container \
  parsingCLIP \
  --db-file /data/gencode.v48.basic.annotation.db \
  --clip-files /data/A27_13.pool_final.NM.bam \
  --out-file /data/A27_13.pool_final.NM.bam.out \
  --genome-dir /data/hg38.ucsc.fasta.chr

# 4. Running omniCLIP

docker run --rm -it \
  -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/pooled_omniclip:/data" \
  kalasnty/omniclip_container:latest run_omniCLIP \
  --db-file /data/gencode.v48.basic.annotation.db \
  --bg-dat /data/A27_13_In.pool_final.NM.bam.out \
  --clip-dat /data/A27_13.pool_final.NM.bam.out \
  --out-dir /data/A27_13.pool_final.NM.omniclip

 # 5. Parse background and CLIP BAMs (individual replicates)

# Background parsing
docker run --rm -it \
  -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/":/data \
  kalasnty/omniclip_container \
  parsingBG \
  --db-file /data/gencode.v48.basic.annotation.db \
  --bg-files /data/A27_13_1In_2.CLAM.genomicAllAligned.NM.bam \
  --out-file /data/A27_13_1In_2.CLAM.genomicAllAligned.NM.bam.out \
  --genome-dir /data/hg38.ucsc.fasta.chr

docker run --rm -it \
  -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/":/data \
  kalasnty/omniclip_container \
  parsingBG \
  --db-file /data/gencode.v48.basic.annotation.db \
  --bg-files /data/A27_13_2In_2.CLAM.genomicAllAligned.NM.bam \
  --out-file /data/A27_13_2In_2.CLAM.genomicAllAligned.NM.bam.out \
  --genome-dir /data/hg38.ucsc.fasta.chr

# CLIP parsing
docker run --rm -it \
  -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/":/data \
  kalasnty/omniclip_container \
  parsingCLIP \
  --db-file /data/gencode.v48.basic.annotation.db \
  --clip-files /data/A27_13_1_2.CLAM.genomicAllAligned.NM.bam \
  --out-file /data/A27_13_1_2.CLAM.genomicAllAligned.NM.bam.out \
  --genome-dir /data/hg38.ucsc.fasta.chr

docker run --rm -it \
  -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/":/data \
  kalasnty/omniclip_container:latest parsingCLIP \
  --db-file /data/gencode.v48.basic.annotation.db \
  --clip-files /data/A27_13_2_2.CLAM.genomicAllAligned.NM.bam \
  --out-file /data/A27_13_2_2.CLAM.genomicAllAligned.NM.bam.out \
  --genome-dir /data/hg38.ucsc.fasta.chr


# 6. Repair HDF5 files if needed (changes in the file names)

# Replicate 1
h5repack /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_1_2.h5 \
         /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_1_2_repaired.h5

h5repack /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_1In_2_bg.h5 \
         /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_1In_2_bg_repaired.h5

# Replicate 2
h5repack /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_2_2.h5 \
         /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_2_2_repaired.h5

h5repack /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_2In_2_bg.h5 \
         /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_2In_2_bg_repaired.h5

# Validation checks

h5dump -H /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_2_2_repaired.h5
h5dump -H /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/A27_13_2In_2_bg_repaired.h5
ls -lh /home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/*repaired.h5

# 7. Run omniCLIP on individual replicates (changes in the file names)

# Replicate 1 with repaired files
docker run --rm -it \
  -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/":/data \
  kalasnty/omniclip_container:latest run_omniCLIP \
  --db-file /data/gencode.v48.basic.annotation.db \
  --bg-dat /data/A27_13_1In_2_bg_repaired.h5 \
  --clip-dat /data/A27_13_1_2_repaired.h5 \
  --out-dir /data/A27_13_1_2.h5.out

# Replicate 2 with repaired files
docker run --rm -it \
  -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/":/data \
  kalasnty/omniclip_container:latest run_omniCLIP \
  --db-file /data/gencode.v48.basic.annotation.db \
  --bg-dat /data/A27_13_2In_2_bg_repaired.h5 \
  --clip-dat /data/A27_13_2_2_repaired.h5 \
  --out-dir /data/A27_13_2_2.h5.out

# 8. Run OmniCLIP with --ign-diag (ignoring the diagnostic event) : 

docker run --rm -it -v "/home/tanasa/16_hours_star_align_richi/A27_13_star_clam_bam/pooled_omniclip:/data" \
kalasnty/omniclip_container:latest run_omniCLIP \
--db-file /data/gencode.v48.basic.annotation.db \
--bg-dat /data/A27_13_In.pool_final.NM.bam.out \
--clip-dat /data/A27_13.pool_final.NM.bam.out \
--out-dir /data/A27_13.pool_final.NM.omniclip.ign.diag \
--ign-diag \
--skip_diag_event_mdl
