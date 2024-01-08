# Usefult Github, sourceforge, and other sites
https://github.com/voutcn/megahit
https://github.com/ablab/spades
https://github.com/BinPro/CONCOCT
https://sourceforge.net/projects/maxbin2/
https://github.com/bxlab/metaWRAP
https://anvio.org/
https://github.com/AnantharamanLab/VIBRANT

#---------------------------------------------- Quality Control of raw reads --------------------------------------------------------#

  # QC was completed for all raw read files ussing bbmap bbduk the flags in the code are explained below

  # Flags
    # first bbduk.sh command line (lines 31-33)
      # ktrim=rl --> trim to the left and right 
      # k=23 --> Kmer length used for finding contaminants
      # mink=11 --> looks for shorter kmers at read tips down to length specified
      # hdist=1 --> Maximum Hamming distance for ref kmers
      # tpe --> when kmer right trimming, trim both reads to the minimum length together
      # tbo --> trim adapters based on where paired reads overlap
      # qtrim=rl --> trim read ends to remove bases with quality below trimq
      # trimq=30 --> regiond with average quality below speicified value will be trimmed
      # ref=/nethome/cxs1538/src/bbmap/resources/adapters.fa --> refrence fasta file with adapters
    # second bbduk.sh command line (lines 34-35)
      # maq= 30 --> reads with average quality (after trimming) below this value will be discarded
    # third bbduk.sh command line (lines 37-39)
      # ref=/nethome/cxs1538/src/bbmap/resources/phix174_ill.ref.fa.gz --> refrence database for phix
      # k=31 --> Kmer length used for finding contaminants
      # hdist=1 --> Maximum Hamming distance for ref kmers
    # fourth bbduk.sh command line (lines 43-44)
      # entropy=0.90 --> filter out reads that have entropy below this value 

  #!/bin/bash 
  for f in *_R1_001.fastq; do 
          echo ${f} 
          name=$(echo ${f} | sed 's/_R1_001.fastq//g') 
          bsub  echo ${name} 
          /nethome/cxs1538/src/bbmap/bbduk.sh -Xmx512m -da in1=${name}_R1_001.fastq in2=${name}_R2_001.fastq \
          out1=${name}_R1_001_out1.fastq out2=${name}_R2_001_out1.fastq \
          ktrim=rl k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=30 ref=/nethome/cxs1538/src/bbmap/resources/adapters.fa 
          /nethome/cxs1538/src/bbmap/bbduk.sh -Xmx512m -da in1=${name}_R1_001_out1.fastq in2=${name}_R2_001_out1.fastq  \
          out1=${name}_R1_001_out2.fastq out2=${name}_R2_001_out2.fastq maq=30 
          /nethome/cxs1538/src/bbmap/bbduk.sh -Xmx512m -da in1=${name}_R1_001_out2.fastq in2=${name}_R2_001_out2.fastq \
          out1=${name}_R1_001_out3.fastq out2=${name}_R2_001_out3.fastq \
          ref=/nethome/cxs1538/src/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 
          /nethome/cxs1538/src/bbmap/bbduk.sh -Xmx512m -da in1=${name}_R1_001_out3.fastq in2=${name}_R2_001_out3.fastq \
          out1=${name}_R1_001_outfinal.fastq out2=${name}_R2_001_outfinal.fastq entropy=0.90  
  done 

  # the quality of the all the reads was visualized using fastqc ls 
  for f in *fastq; do
  	fastqc $f
  done 
  # the seawater metagenome and virome reads (forward and reverse) were good to continue with
  # the sargassum metagenome forward reads were also good to continue with

  # the sargassum metagenome reverse reads had at 15bp region that was a lot of G/A repeats so I
  # removed the first 15 bp of each sargassum reverse read

  #!/bin/bash
  eval "$(conda shell.bash hook)"
  for f in *_R2_001_outfinal.fastq; do
          sample=$(echo ${f} | sed 's/_R2_001_outfinal.fastq//g')
          fastp  -I "$sample"_R2_001_outfinal.fastq -O "$sample"_R2_trim.fastq -f 15 --disable_adapter_trimming
  done

#------------------------------------------ Read taxonomic profiles for metagenomes with Kaiju ----------------------------------------------------#
  # Sargassum
  #!/bin/bash
  /projectnb/viralecology/databases/kaiju/src/kaiju-multi -z 25 \
  -t /projectnb/viralecology/databases/kaiju-db/progenomes_100922/nodes.dmp \
  -f /projectnb/viralecology/databases/kaiju-db/progenomes_100922/kaiju_db_progenomes.fmi \
  -i Sargasso-1_R1_001_outfinal.fastq,Sargasso-2_R1_001_outfinal.fastq,Sargasso-3_R1_001_outfinal.fastq,Sargasso-4_R1_001_outfinal.fastq,Sargasso-6_R1_001_outfinal.fastq \
  -j Sargasso-1_R2_trim_outfinal.fastq,Sargasso-2_R2_trim_outfinal.fastq,Sargasso-3_R2_trim_outfinal.fastq,Sargasso-4_R2_trim_outfinal.fastq,Sargasso-6_R2_trim_outfinal.fastq \
  -o Sargasso-1.out,Sargasso-2.out,Sargasso-3.out,Sargasso-4.out,Sargasso-6.out

  #!/bin/bash
  for f in *.out ;do
          sample=$(echo ${f} | sed 's/.out//g')
          /projectnb/viralecology/databases/kaiju/src/kaiju2table -t /projectnb/viralecology/databases/kaiju-db/progenomes_100922/nodes.dmp \
          -n /projectnb/viralecology/databases/kaiju-db/progenomes_100922/names.dmp -u -r order -l phylum,class,order \
          -o "$sample"_phylum_class.tsv $sample.out
  done

  # Seawater Cellular Fraction
  #!/bin/bash
  /projectnb/viralecology/databases/kaiju/src/kaiju-multi -z 25 \
  -t /projectnb/viralecology/databases/kaiju-db/progenomes_100922/nodes.dmp \
  -f /projectnb/viralecology/databases/kaiju-db/progenomes_100922/kaiju_db_progenomes.fmi \
  -i Seawater-Sterivex-1.R1_001_outfinal.fastq,Seawater-Sterivex-2.R1_001_outfinal.fastq,Seawater-Sterivex-3.R1_001_outfinal.fastq,Seawater-Sterivex-4.R1_001_outfinal.fastq,Seawater-Sterivex-5.R1_001_outfinal.fastq \
  -j Seawater-Sterivex-1.R2_001_outfinal.fastq,Seawater-Sterivex-2.R2_001_outfinal.fastq,Seawater-Sterivex-3.R2_001_outfinal.fastq,Seawater-Sterivex-4.R2_001_outfinal.fastq,Seawater-Sterivex-5.R2_001_outfinal.fastq  \
  -o Seawater-Sterivex-1.out,Seawater-Sterivex-2.out,Seawater-Sterivex-3.out,Seawater-Sterivex-4.out,Seawater-Sterivex-5.out

  #!/bin/bash
  for f in *.out ;do
          sample=$(echo ${f} | sed 's/.out//g')
          /projectnb/viralecology/databases/kaiju/src/kaiju2table -t /projectnb/viralecology/databases/kaiju-db/progenomes_100922/nodes.dmp \
          -n /projectnb/viralecology/databases/kaiju-db/progenomes_100922/names.dmp -u -r order -l order \
          -o "$sample"_order.tsv $sample.out
  done 


#-------------------------------------------------------- Contig Assembly ----------------------------------------------------------------------#
  # Sargassum Co-assemblies
  nethome/aks206/miniconda3/envs/utils/bin/megahit --presets meta-large \
  -1 Sargasso-1_R1_001_outfinal.fastq,Sargasso-2_R1_001_outfinal.fastq,Sargasso-3_R1_001_outfinal.fastq,Sargasso-4_R1_001_outfinal.fastq,Sargasso-6_R1_001_outfinal.fastq \
  -2 Sargasso-1_R2_trim.fastq,Sargasso-2_R2_trim.fastq,Sargasso-3_R2_trim.fastq,Seawater-Sargasso-4_R2_trim.fastq,Sargasso-6_R2_trim.fastq \
  --min-contig-len 1000 -m 0.85 \
  -o Sargasso/megahit/ -t 16

  # Seawater Cellular Fraction Co-assemblies 
  /nethome/aks206/miniconda3/envs/utils/bin/megahit --presets meta-large \
  -1 Seawater-Sterivex-1_R1_001_outfinal.fastq,Seawater-Sterivex-2_R1_001_outfinal.fastq,Seawater-Sterivex-3_R1_001_outfinal.fastq,Seawater-Sterivex-4_R1_001_outfinal.fastq,Seawater-Sterivex-5_R1_001_outfinal.fastq \
  -2 Seawater-Sterivex-1_R2_001_outfinal.fastq,Seawater-Sterivex-2_R2_001_outfinal.fastq,Seawater-Sterivex-3_R2_001_outfinal.fastq,Seawater-Sterivex-4_R2_001_outfinal.fastq,Seawater-Sterivex-5_R2_001_outfinal.fastq \
  --min-contig-len 1000 -m 0.85 \
  -o Seawater-Sterivex/megahit/ -t 16

# Per Sample Assemby for Seawater Viral Fraction
  for f in *R1_001_outfinal.fastq; do 
    name=$(echo ${f} | sed 's/R1_001_outfinal.fastq//g') 
    mkdir ${name}_metaspades 
    /nethome/aks206/mybin/SPAdes-3.15.4-Linux/bin/spades.py --meta --pe1-1 ${name}R1_001_outfinal.fastq \
    --pe1-2 ${name}R2_001_outfinal.fastq \
    --only-assembler -o ${name}_metaspades
  done 
  # concatenate everything into one file 
  for f in *_metaspades ;do
    name=$(echo ${f} | sed 's/_metaspades//g')
    cd $f
    mv contigs.fasta ${name}_contigs.fa
    cp ${name}_contigs.fa /scratch/projects/sarg2022/Seawater-virome/spades 
  done
  cat *_contigs.fa > all_Seawater-virome_contigs_metaspades.fasta

#----------------------------------------- bMAG Assembly with CONCOCT, MaxBin2, and MetaBAT2 -----------------------------------------------------------#
  # mapping reads to contigs for coverage data

  # Sargassum
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/anvio-7.1
  bowtie2-build /scratch/projects/sarg2022/Sargasso/megahit/final.contigs.fa Sargasso_bowtieDB
  conda deactivate

  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/anvio-7.1
  for f in *_R1_001_outfinal.fastq; do
    sample=$(echo ${f} | sed 's/_R1_001_outfinal.fastq//g')
    bowtie2 --threads 20 -x /scratch/projects/sarg2022/Sargasso/mapping/Sargasso_bowtieDB -1 "$sample"_R1_001_outfinal.fastq -2 "$sample"_R2_trim.fastq \
    -S /scratch/projects/sarg2022/Sargasso/mapping/"$sample".sam
    samtools view -S -b /scratch/projects/sarg2022/Sargasso/mapping/"$sample".sam | samtools sort > /scratch/projects/sarg2022/Sargasso/mapping/"$sample".sorted.bam
    samtools index /scratch/projects/sarg2022/Sargasso/mapping/"$sample".sorted.bam
  done
  conda deactivate
  conda activate /nethome/pxh390/miniconda3/envs/metawrap-env
  cd /scratch/projects/sarg2022/Sargasso/mapping/
  for f in * .sam; do
    sample=$(echo ${f} | sed 's/.sam//g')
    jgi_summarize_bam_contig_depths --outputDepth "$sample".alignment.sorted.depth.txt $samp.sorted.bam
    pileup.sh in="$sample".sorted.bam out="$sample".alignment.cov.txt overwrite=true
    awk '{print $1"\t"$5}' "$sample".alignment.cov.txt | grep -v '^#' > "$sample".alignment.abundance.txt
  done
  conda deactivate 

  # Seawater Cellular Fraction
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/anvio-7.1
  bowtie2-build /scratch/projects/sarg2022/Seawater-Sterivex/megahit/final.contigs.fa Seawater-Sterivex_bowtieDB
  conda deactivate

  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/anvio-7.1
  for f in *_R1_001_outfinal.fastq; do
    sample=$(echo ${f} | sed 's/_R1_001_outfinal.fastq//g')
    bowtie2 --threads 20 -x /scratch/projects/sarg2022/Sargasso/mapping/Seawater-Sterivex_bowtieDB -1 "$sample"_R1_001_outfinal.fastq -2 "$sample"_R2_trim.fastq \
    -S /scratch/projects/sarg2022/Seawater-Sterivex/mapping/"$sample".sam
    samtools view -S -b /scratch/projects/sarg2022/Seawter-Sterivex/mapping/"$sample".sam | samtools sort > /scratch/projects/sarg2022/Seawater-Sterivex/mapping/"$sample".sorted.bam
    samtools index /scratch/projects/sarg2022/Sargasso/Seawater-Sterivex/"$sample".sorted.bam
  done
  conda deactivate
  conda activate /nethome/pxh390/miniconda3/envs/metawrap-env
  cd /scratch/projects/sarg2022/Seawater-Sterivex/mapping/
  for f in * .sam; do
    sample=$(echo ${f} | sed 's/.sam//g')
    jgi_summarize_bam_contig_depths --outputDepth "$sample".alignment.sorted.depth.txt $samp.sorted.bam
    pileup.sh in="$sample".sorted.bam out="$sample".alignment.cov.txt overwrite=true
    awk '{print $1"\t"$5}' "$sample".alignment.cov.txt | grep -v '^#' > "$sample".alignment.abundance.txt
  done
  conda deactivate

  # Binning with CONCOCT

  # Sargassum
  cd /scratch/projects/sarg2022/Sargasso
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate ~/miniconda3/envs/metawrap-env
    cut_up_fasta.py megahit/final.contigs.fa -c 10000 -o 0 \
                  --merge_last -b binning/contigs_10K.bed > binning/contigs_10K.fa
    #
    ls mapping/*.bam | sed 's/.sorted.bam//g' | sed 's/mapping\///g' > samp-names.tsv
    BAMS=`ls mapping/*.bam | python -c 'import sys; print(" ".join([x.strip() for x in sys.stdin.readlines()]))'`
    concoct_coverage_table.py binning/contigs_10K.bed --samplenames samp-names.tsv \
        $BAMS > binning/coverage_table.csv
    # run concoct
    concoct --composition_file binning/contigs_10K.fa \
            --coverage_file binning/coverage_table.csv -b binning/
    # Merge subcontig clustering into original contig clustering
    merge_cutup_clustering.py binning/clustering_gt1000.csv > binning/clustering_merged.csv
    # Extract bins as individual FASTA
    extract_fasta_bins.py megahit/final.contigs.fa binning/clustering_merged.csv \
            --output_path binning/concoct/
            # CONCOCT output the bin name as a "number".fa (e.g. "1.fa")
            # which a lot of software is NOT going to like filenames starting
            # number -_- (change to "bin1.fa")
            for f in binning/concoct/*.fa; do
              mv $f binning/concoct/bin"$(basename $f)"
            done
  # Seawater Cellular Fraction
  cd /scratch/projects/sarg2022/Seawater-Sterivex
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate ~/miniconda3/envs/metawrap-env
    cut_up_fasta.py megahit/final.contigs.fa -c 10000 -o 0 \
                  --merge_last -b binning/contigs_10K.bed > binning/contigs_10K.fa
    #
    ls mapping/*.bam | sed 's/.sorted.bam//g' | sed 's/mapping\///g' > samp-names.tsv
    BAMS=`ls mapping/*.bam | python -c 'import sys; print(" ".join([x.strip() for x in sys.stdin.readlines()]))'`
    concoct_coverage_table.py binning/contigs_10K.bed --samplenames samp-names.tsv \
        $BAMS > binning/coverage_table.csv
    # run concoct
    concoct --composition_file binning/contigs_10K.fa \
            --coverage_file binning/coverage_table.csv -b binning/
    # Merge subcontig clustering into original contig clustering
    merge_cutup_clustering.py binning/clustering_gt1000.csv > binning/clustering_merged.csv
    # Extract bins as individual FASTA
    extract_fasta_bins.py megahit/final.contigs.fa binning/clustering_merged.csv \
            --output_path binning/concoct/
            # CONCOCT output the bin name as a "number".fa (e.g. "1.fa")
            # which a lot of software is NOT going to like filenames starting
            # number -_- (change to "bin1.fa")
            for f in binning/concoct/*.fa; do
              mv $f binning/concoct/bin"$(basename $f)"
            done
  conda deactivate
  # Binning with MetaBAT2

    # Sargassum
    cd /scratch/projects/sarg2022/Sargasso
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate ~/miniconda3/envs/metawrap-env
    jgi_summarize_bam_contig_depths --outputDepth mapping/depth.txt mapping/*.bam
    metabat2 -m 1500 -t 8 -i megahit/final.contigs.fa \
            -a mapping/depth.txt -o binning/metabat2/bin
    conda deactivate 

    # Seawater Cellular Fraction
    cd /scratch/projects/sarg2022/Seawater-Sterivex
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate ~/miniconda3/envs/metawrap-env
    jgi_summarize_bam_contig_depths --outputDepth mapping/depth.txt mapping/*.bam
    metabat2 -m 1500 -t 8 -i megahit/final.contigs.fa \
            -a mapping/depth.txt -o binning/metabat2/bin
    conda deactivate 

  # Binning with MaxBIN2

    # Sargassum
      cd /scratch/projects/sarg2022/Sargasso
      ls mapping/*.alignment.abundance.txt > mapping/abund.list.txt
      run_MaxBin.pl -thread 4 -contig megahit/final.contigs.fa \
                    -out binning/maxbin2/bin \
                    -abund_list mapping/abund.list.txt
      # mega annoying things about maxbin is that it outputs bins as .fasta
      # while the other binners output to .fa --> MetaWRAP DOESNT like this...
      # so need to change all the names
      mv binning/maxbin2/*.fasta binning/maxbin2/bin
      for f in binning/maxbin2/bin/*.fasta; do
      samp=$(basename $f .fasta).fa
      mv $f binning/maxbin2/bin/"$samp"
      done

    # Seawater Cellular Fraction 
      cd /scratch/projects/sarg2022/Seawater-Sterivex
      ls mapping/*.alignment.abundance.txt > mapping/abund.list.txt
      run_MaxBin.pl -thread 4 -contig megahit/final.contigs.fa \
                    -out binning/maxbin2/bin \
                    -abund_list mapping/abund.list.txt
      # mega annoying things about maxbin is that it outputs bins as .fasta
      # while the other binners output to .fa --> MetaWRAP DOESNT like this...
      # so need to change all the names
      mv binning/maxbin2/*.fasta binning/maxbin2/bin
      for f in binning/maxbin2/bin/*.fasta; do
      samp=$(basename $f .fasta).fa
      mv $f binning/maxbin2/bin/"$samp"
      done
#----------------------------------------- bMAG refinement, reassembly, taxonomy, and quantification with MetaWrap -----------------------------------------------------------#
  # Bin Refinement
  # Sargassum 
  cd /scratch/projects/sarg2022/Sargasso
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/metawrap-env
  metawrap bin_refinement -o binning/refined_bins_20-10 -c 20 -x 10 -t 4 \
                          -A binning/metabat2/ \
                          -B binning/maxbin2/bin \
                          -C binning/concoct/
  conda deactivate

  # Seawater Cellular Fraction 
  cd /scratch/projects/sarg2022/Seawater-Sterivex
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/metawrap-env
  metawrap bin_refinement -o binning/refined_bins_20-10 -c 20 -x 10 -t 4 \
                          -A binning/metabat2/ \
                          -B binning/maxbin2/bin \
                          -C binning/concoct/
  conda deactivate

  # Bin Re-assembly
  # Sargasum
  cd /scratch/projects/sarg2022/Sargasso
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/metawrap-env
  metawrap reassemble_bins -o /scratch/projects/sarg2022/Sargasso/BIN_REASSEMBLY \
  -1 fastq/ALL_READS_1.fastq -2 fastq/ALL_READS_2.fastq \
  -t 96 -m 800 -c 50 -x 10 -b binning/refined_bins_20-10/metawrap_20_10_bins/
  conda deactivate

  # Seawater Cellular Fraction
  cd /scratch/projects/sarg2022/Seawater-Sterivex
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/metawrap-env
  metawrap reassemble_bins -o /scratch/projects/sarg2022/Seawater-Sterivex/BIN_REASSEMBLY \
  -1 fastq/ALL_READS_1.fastq -2 fastq/ALL_READS_2.fastq \
  -t 96 -m 800 -c 50 -x 10 -b binning/refined_bins_20-10/metawrap_20_10_bins/
  conda deactivate 

  # Bin taxonomy
  # Sargassum
  cd /scratch/projects/sarg2022/Sargasso
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/metawrap-env
  metawrap classify_bins -b /scratch/projects/sarg2022/Sargasso/BIN_REASSEMBLY -o BIN_CLASSIFICATION -t 48
  conda deactivate

  # Seawwater Cellular Fraction 
  cd /scratch/projects/sarg2022/Seawater-Sterivex
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/metawrap-env
  metawrap classify_bins -b scratch/projects/sarg2022/Seawater-Sterivex/BIN_REASSEMBLY -o BIN_CLASSIFICATION -t 48
  conda deactivate

  # Quantification of NON-REASSEMBLED bins with Metawrap (metawrap suggests not quantifying reassembled bins)
    # Sargassum
    cd /scratch/projects/sarg2022/Sargasso
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/aks206/miniconda3/envs/metawrap-env
    metawrap quant_bins -b NON_REASSEMBLED_Dereplicated_50_10_MAGs -o NON_REASSEMBLED_50_10_Quant \
     -a megahit/final_contigs_fixed.fa fastq/metawrap_reads/*.fastq
    conda deactivate

    # Seawater Cellular Fraction 
    cd scratch/projects/sarg2022/Seawater-Sterivex
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/aks206/miniconda3/envs/metawrap-env
    metawrap quant_bins -b NON_REASSEMBLED_Dereplicated_50_10_MAGs -o NON_REASSEMBLED_50_10_MAG_Quant \
     -a megahit/final_contigs_fixed.fa fastq/metawrap_reads/*.fastq
    conda deactivate

#--------------------------------------- Generation of bMAG Pecernt Nucleotide Identity Tree ----------------------------------#

  # Dereplication of >=50% complete and =<10% contamination bins with Anvio and fastANI
    # Sargassum
    # create an anvio contig database for each bin
    cd Sargasso/BIN_REASSEMBLY/50_10_bins
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/aks206/miniconda3/envs/anvio-7.1
    for f in *.fa; do
            bin=$(echo ${f} | sed 's/.fa//g')
            anvi-gen-contigs-database -T 16 -f $f \
            -o Sargasso/bin_dereplication/"$bin"_contig.db
    done
    conda deactivate

    # make the external genome file
    cd Sargasso/BIN_REASSEMBLY/50_10_bins
    # get bin names
    ls | sed -e 's/.fa//g' > bin_names.txt
    sed -e 's/bin_names.txt//g' | sed '/^$/d' > bin_names1.txt
    cd Sargasso/bin_dereplication
    # get contig db paths
    realpath *.db > paths.txt
    # create column header
    echo 'name' > name.txt
    echo 'contigs_db_path' > contig_path.txt
    # paste everything togeth
    paste name.txt contig_path.txt > name_contigs_path.txt
    paste bin_names1.txt paths.txt > bin_names_paths.txt
    cat name_contigs_path bin_names_path.txt > external-genome.txt

    # run fastANI through Anvio
    cd Sargasso/bin_dereplication
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/pxh390/miniconda3/envs/anvio-7.1/
    anvi-dereplicate-genomes -e external-genome.txt \
            --program fastANI \
            -o fastANI_derep95PERC \
            --min-full-percent-identity 0.95 \
            --similarity-threshold 0.95 \
            --min-fraction 25
    conda deactivate

    # Seawater Cellular Fraction
    # create an anvio contig database for each bin
    cd Seawater-Sterivex/BIN_REASSEMBLY/50_10_bins
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/aks206/miniconda3/envs/anvio-7.1
    for f in *.fa; do
            bin=$(echo ${f} | sed 's/.fa//g')
            anvi-gen-contigs-database -T 16 -f $f \
            -o Seawater-Sterivex/bin_dereplication/"$bin"_contig.db
    done
    conda deactivate

    # make the external genome file
    cd Seawater-Sterivex/BIN_REASSEMBLY/50_10_bins
    # get bin names
    ls | sed -e 's/.fa//g' > bin_names.txt
    sed -e 's/bin_names.txt//g' | sed '/^$/d' > bin_names1.txt
    cd Seawater-Sterivex/bin_dereplication
    # get contig db paths
    realpath *.db > paths.txt
    # create column header
    echo 'name' > name.txt
    echo 'contigs_db_path' > contig_path.txt
    # paste everything togeth
    paste name.txt contig_path.txt > name_contigs_path.txt
    paste bin_names1.txt paths.txt > bin_names_paths.txt
    cat name_contigs_path bin_names_path.txt > external-genome.txt

    # run fastANI through Anvio
    cd Seawater-Sterivex/bin_dereplication
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/pxh390/miniconda3/envs/anvio-7.1/
    anvi-dereplicate-genomes -e external-genome.txt \
            --program fastANI \
            -o fastANI_derep95PERC \
            --min-full-percent-identity 0.95 \
            --similarity-threshold 0.95 \
            --min-fraction 25
    conda deactivate

#--------------------------------------- Viral contigs indetification and prophage confirmation with VIBRANT and CheckV ----------------------------------#
  # Identification of viral contigs and bacterial flanking regions in bMAGs that >=50% complete and =<10% contamination with VIBRANT
    # Sargassum
    cd Sargasso/bin_dereplication/fastANI_derep95PERC/GENOMES
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/aks206/miniconda3/envs/VIBRANT
    for f in *.fa ;do
      bin=$(echo ${f} | sed 's/.fa//g')
      /nethome/aks206/miniconda3/envs/VIBRANT/VIBRANT/VIBRANT_run.py -i $bin.fa -o bin_VIBRANT/"$bin"_VIBRANT
    done 
    conda deactivate

    # Seawater Cellular Fraction
    cd Seawater-Sterivex/bin_dereplication/fastANI_derep95PERC/GENOMES
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/aks206/miniconda3/envs/VIBRANT
    for f in *.fa ;do
      bin=$(echo ${f} | sed 's/.fa//g')
      /nethome/aks206/miniconda3/envs/VIBRANT/VIBRANT/VIBRANT_run.py -i $bin.fa -o bin_VIBRANT/"$bin"_VIBRANT
    done 
    conda deactivate

  # Identification of bacterial flanking regions with CheckV
    # Sargassum
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/nsv19/anaconda3/envs/vRhyme
    checkv=/nethome/nsv19/anaconda3/envs/vRhyme/bin/checkv
    checkv end_to_end /scratch/projects/sarg2022/Sargasso/megahit/final.contigs.fa CheckV \
    -d /nethome/nsv19/anaconda3/envs/vRhyme/checkv-db-v1.4 -t 16
    conda deactivate

    #Seawater Cellular Fraction
    #!/bin/bash
    eval "$(conda shell.bash hook)"
    conda activate /nethome/nsv19/anaconda3/envs/vRhyme
    checkv=/nethome/nsv19/anaconda3/envs/vRhyme/bin/checkv
    checkv end_to_end /scratch/projects/sarg2022/Seawater-Sterivex/megahit/final.contigs.fa CheckV \
    -d /nethome/nsv19/anaconda3/envs/vRhyme/checkv-db-v1.4 -t 16
    conda deactivate

#--------------------------------------- Identification of viral fragments in all contigs with VIBRANT ----------------------------------#
  # dereplicate contigs with cdhit at 98 ANI
  # Sargassum
  /nethome/ask206/mybin/cd-hit-v4.8.1-2019-0228/cd-hit-est -i Sargasso/megahit/final_contigs.fa \
  -o Sargassum_contigs_cdhit98.fa -n 10 -M 20000 -c 0.98 
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/VIBRANT
  /nethome/aks206/mybin/VIBRANT/VIBRANT_run.py -i Sargasso/megahit/Sargassum_contigs_cdhit98.fa \
  -o Sargasso/Sargassum_contigs_cdhit_98_VIBRANT
  
  # Seawater Cellular Fraction
  /nethome/ask206/mybin/cd-hit-v4.8.1-2019-0228/cd-hit-est -i Seawater-Sterivex/megahit/final_contigs.fa \
  -o Seawater-Sterivex_contigs_cdhit98.fa -n 10 -M 20000 -c 0.98 
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/VIBRANT
  /nethome/aks206/mybin/VIBRANT_run.py -i Seawater-Sterivex/megahit/Seawater-Sterivex_contigs_cdhit98.fa \
  -o Seawater-Sterivex/Seawater-Sterivex_contigs_cdhit_98_VIBRANT
  
  # Seawater-virome
  /nethome/ask206/mybin/cd-hit-v4.8.1-2019-0228/cd-hit-est -i Seawater-virome/megahit/final_contigs.fa \
  -o Sargassum_contigs_cdhit98.fa -n 10 -M 20000 -c 0.98
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/aks206/miniconda3/envs/VIBRANT

#------------------------------------------------- vMAG assembly with vRhyme ----------------------------------------------#
  # Sargassum
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/nsv19/anaconda3/envs/vRhyme
  vRhyme -i /scratch/projects/sarg2022/Sargasso/VIBRANT_Sargasso_contigs_cdhit/VIBRANT_phages_Sargasso_contigs_cdhit/Sargasso_contigs_cdhit.phages_combined.fna \
  -o /scratch/projects/sarg2022/Sargasso/Viral_binning/All_Sargassum_vRhyme_output \
  -r /scratch/projects/sarg2022/Sargasso/fastq/ALL_READS* \
  --method longest
  conda deactivate

  #Seawater Cellular Fraction 
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/nsv19/anaconda3/envs/vRhyme
  vRhyme -i /scratch/projects/sarg2022/Seawater-Sterivex/VIBRANT_Seawater-Sterivex_contigs_cdhit/VIBRANT_phages_Seawater-Sterivex_contigs_cdhit/Seawater-Sterivex_contigs_cdhit.phages_combined.fna \
  -o /scratch/projects/sarg2022/Seawater-Sterivex/Viral_binning/All_Seawater-Sterivex_vRhyme_output \
  -r /scratch/projects/sarg2022/Seawater-Sterivex/fastq/ALL_READS* \
  --method longest
  conda deactivate 

  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/nsv19/anaconda3/envs/vRhyme
  vRhyme -i /scratch/projects/sarg2022/Seawater-virome/VIBRANT_all_sw_virome_contigs_fixed/VIBRANT_phages_all_sw_virome_contigs_fixed/Seawater-virome_contigs_fixed.phages_combined.fna \
  -o /scratch/projects/sarg2022/Seawater-virome/Viral_binning/All_Seawater-virome_vRhyme_output \
  -r /scratch/projects/sarg2022/Seawater-virome/fastq/ALL_READS* \
  --method longest
  conda deactivate

#------------------------------------------------- vMAG Dereplication with virathon ---------------------------------------------------#
  # Remove duplicate contigs in vMAG

  # Run a blast to identify identical contigs
  #Sargassum
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate utils
  for f in *.fasta; do
          sname=$(echo ${f} | sed 's/.fasta//g')
    makeblastdb -in ${f} -dbtype nucl -out /scratch/projects/sarg2022/vMAG_dereplication/blastn_pre-derep_db/Sargassum_dbs/${sname}_db
    blastn -db /scratch/projects/sarg2022/vMAG_dereplication/blastn_pre-derep_db/Sargassum_dbs/${sname}_db \
    -query /scratch/projects/sarg2022/vMAG_dereplication/Sargassum_vMAGs/${f} \
    -out /scratch/projects/sarg2022/vMAG_dereplication/blastn_out_97id95cov/Sargassum/${sname}.out \
    -outfmt "6 qseqid sseqid pident qcovs qlen slen length mismatch qstart qend sstart send evalue" -num_threads 16 
  done
  conda deactivate 

  #Seawater Cellular Fraction
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate utils
  for f in *.fasta; do
          sname=$(echo ${f} | sed 's/.fasta//g')
    makeblastdb -in ${f} -dbtype nucl -out /scratch/projects/sarg2022/vMAG_dereplication/blastn_pre-derep_db/Seawater-Sterivex_dbs/${sname}_db
    blastn -db /scratch/projects/sarg2022/vMAG_dereplication/blastn_pre-derep_db/Seawater-Sterivex_dbs/${sname}_db \
    -query /scratch/projects/sarg2022/vMAG_dereplication/Seawater-Sterivex_vMAGs/${f} \
    -out /scratch/projects/sarg2022/vMAG_dereplication/blastn_out_97id95cov/Seawater-Sterivex/${sname}.out \
    -outfmt "6 qseqid sseqid pident qcovs qlen slen length mismatch qstart qend sstart send evalue" -num_threads 16 
  done
  conda deactivate 

  # Seawater viral fraction
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate utils
  for f in *.fasta; do
          sname=$(echo ${f} | sed 's/.fasta//g')
    makeblastdb -in ${f} -dbtype nucl -out /scratch/projects/sarg2022/vMAG_dereplication/blastn_pre-derep_db/Seawater-virome_dbs/${sname}_db
    blastn -db /scratch/projects/sarg2022/vMAG_dereplication/blastn_pre-derep_db/Seawater-virome_dbs/${sname}_db \
    -query /scratch/projects/sarg2022/vMAG_dereplication/Seawater-virome_vMAGs/${f} \
    -out /scratch/projects/sarg2022/vMAG_dereplication/blastn_out_97id95cov/Seawater-virome/${sname}.out \
    -outfmt "6 qseqid sseqid pident qcovs qlen slen length mismatch qstart qend sstart send evalue" -num_threads 16 
  done
  conda deactivate 

  # Get the list of sequence IDs for the contigs that are matches (>97 percent identity, >95% coverage)
  # And retain the longest contig
  # Sargassum
  for f in *.out; do
    name=$(echo ${f} | sed 's/.out//g')
    awk '{ if ($3>97) print }' ${name}.out \
    | awk '{ if ($1!=$2) print}' | awk '{ if ($5<$6) print}' |  awk '{ if ($4>95) print }' | awk '{print $1}' | sort | uniq \
    > list_to_remove/${name}.remove.list
  done

  # Seawater Cellular Fraction
  for f in *.out; do
    name=$(echo ${f} | sed 's/.out//g')
    awk '{ if ($3>97) print }' ${name}.out \
    | awk '{ if ($1!=$2) print}' | awk '{ if ($5<$6) print}' |  awk '{ if ($4>95) print }' | awk '{print $1}' | sort | uniq \
    > list_to_remove/${name}.remove.list
  done

  # Seawater Viral Fraction 
  for f in *.out; do
    name=$(echo ${f} | sed 's/.out//g')
    awk '{ if ($3>97) print }' ${name}.out \
    | awk '{ if ($1!=$2) print}' | awk '{ if ($5<$6) print}' |  awk '{ if ($4>95) print }' | awk '{print $1}' | sort | uniq \
    > list_to_remove/${name}.remove.list
  done

  # remove the duplicated contigs with seqkit 
  #Sargassum
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate utils
  for f in *.fasta; do
    name=$(echo ${f} | sed 's/.fasta//g')
    echo $sname
    seqkit grep -v -f /scratch/projects/sarg2022/vMAG_dereplication/blastn_out_97id95cov/Sargassum/list_to_remove/${name}.remove.list ${f} \
    > /scratch/projects/sarg2022/vMAG_dereplication/Sargassum_vMAGs/bins_without_dubs/${name}.no_dups.fasta
  done
  conda deactivate 

  # Seawater Cellular Fraction
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate utils
  for f in *.fasta; do
    name=$(echo ${f} | sed 's/.fasta//g')
    echo $sname
    seqkit grep -v -f /scratch/projects/sarg2022/vMAG_dereplication/blastn_out_97id95cov/Seawater-Sterivex/list_to_remove/${name}.remove.list ${f} \
    > /scratch/projects/sarg2022/vMAG_dereplication/Seawater-Sterivex_vMAGs/bins_without_dubs/${name}.no_dups.fasta
  done
  conda deactivate 

  # Seawater Viral Fraction
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate utils
  for f in *.fasta; do
    name=$(echo ${f} | sed 's/.fasta//g')
    echo $sname
    seqkit grep -v -f /scratch/projects/sarg2022/vMAG_dereplication/blastn_out_97id95cov/Seawater-virome/list_to_remove/${name}.remove.list ${f} \
    > /scratch/projects/sarg2022/vMAG_dereplication/Seawater-virome_vMAGs/bins_without_dubs/${name}.no_dups.fasta
  done
  conda deactivate 

  # Convert files back to single line
  # Sargassum
  #!/bin/bash
  for f in *.fasta ; do
          name=$(echo ${f} | sed 's/.no_dups.fasta//g')
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ${f}\
    > /scratch/projects/sarg2022/vMAG_dereplication/Sargassum_vMAGs/bins_without_dubs_single_line/${name}single_line.fasta
  done

  # Seawater Cellular Fraction
  #!/bin/bash
  for f in *.fasta ; do
          name=$(echo ${f} | sed 's/.no_dups.fasta//g')
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ${f}\
    > /scratch/projects/sarg2022/vMAG_dereplication/Seawater-Sterivex_vMAGs/bins_without_dubs_single_line/${name}single_line.fasta
  done

  # Seawater Viral Fraction
  #!/bin/bash
  for f in *.fasta ; do
          name=$(echo ${f} | sed 's/.no_dups.fasta//g')
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ${f}\
    > /scratch/projects/sarg2022/vMAG_dereplication/Seawater-virome_vMAGs/bins_without_dubs_single_line/${name}single_line.fasta
  done

  # N-link the viral bins without duplicate contigs
  # Sargassum
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/nsv19/anaconda3/envs/vRhyme
  link_bin_sequences.py -i /scratch/projects/sarg2022/vMAG_dereplication/Sargassum_vMAGs/bins_without_dubs_single_line/ \
  -o /scratch/projects/sarg2022/vMAG_dereplication/Sargassum_vMAGs/N_linked_bins_without_dubs_02-08-23 \
  -e fasta \
  -n 1000 \
  -c N 
  conda deactivate

  # Seawater Cellular Fraction
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/nsv19/anaconda3/envs/vRhyme
  link_bin_sequences.py -i /scratch/projects/sarg2022/vMAG_dereplication/Seawater-Sterivex_vMAGs/bins_without_dubs_single_line/ \
  -o /scratch/projects/sarg2022/vMAG_dereplication/Seawater-Sterivex_vMAGs/N_linked_bins_without_dubs_02-08-23 \
  -e fasta \
  -n 1000 \
  -c N 
  conda deactivate 

  # Seawater Viral Fraction
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/nsv19/anaconda3/envs/vRhyme
  link_bin_sequences.py -i /scratch/projects/sarg2022/vMAG_dereplication/Seawater-_vMAGs/bins_without_dubs_single_line/ \
  -o /scratch/projects/sarg2022/vMAG_dereplication/Seawater-virome_vMAGs/N_linked_bins_without_dubs_02-08-23 \
  -e fasta \
  -n 1000 \
  -c N 
  conda deactivate 

  # Combine all the N-linked vMAGS and the non-binned contigs that have medium/high quality 
  # based on genome annotations from VIBRANT into one file and dereplicate across samples
  # The resulting viral genomes from this will make up the SSVdb
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/nsv19/anaconda3/envs/virathon
  python3 /nethome/nsv19/anaconda3/envs/virathon/share/virathon/Virathon.py --genome_files ALL_vMAGS_10000_030923.fasta --make_pops True --threads 24
  conda deactivate

#------------------------------------------------- Fractional Abundance of the SSVdb using SMALT ---------------------------------------------------#
  # create SMALT database for the dereplicated SSVdb
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/pxh390/miniconda3/envs/samtools
  smalt index all_phages_050523 /scratch/projects/sarg2022/vMAG_dereplication/ALL_DEREPLICATED_vMAGs_great10000_050523.fasta
  conda deactivate

  # Map the SSVdb to the Sargassum reads 
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/pxh390/miniconda3/envs/samtools
  for f in *_R1_001_outfinal.fastq; do
    sample=$(echo ${f} | sed 's/_R1_001_outfinal.fastq//g')
    smalt map -n 40 -y 0.95 -o /scratch/projects/sarg2022/Sargasso/phage_abundances/$sample.aln.sam \
    /scratch/projects/sarg2022/vMAG_dereplication/all_phages_050523 \
    "$sample"_R1_001_outfinal.fastq "$sample"_R2_trim_outfinal.fastq
    samtools sort -o /scratch/projects/sarg2022/Sargasso/phage_abundances/$sample.aln.sort.bam /scratch/projects/sarg2022/Sargasso/phage_abundances/$sample.aln.sam
    rm /scratch/projects/sarg2022/Sargasso/phage_abundances/$sample.aln.sam
    samtools index /scratch/projects/sarg2022/Sargasso/phage_abundances/$sample.aln.sort.bam
    samtools idxstats /scratch/projects/sarg2022/Sargasso/phage_abundances/$sample.aln.sort.bam > /scratch/projects/sarg2022/Sargasso/phage_abundances/$sample.idxstats
  done

  # Map the SSVdb to the Seawater Cellular Fraction reads 
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/pxh390/miniconda3/envs/samtools
  for f in *.qssghR1_001_outfinal.fastq; do
    sample=$(echo ${f} | sed 's/.R1_001_outfinal.fastq//g')
    smalt map -n 40 -y 0.95 -o /scratch/projects/sarg2022/Seawater-Sterivex/phage_abundances/$sample.aln.sam \
    /scratch/projects/sarg2022/vMAG_dereplication/all_phages_050523 \
    "$sample"_R1_001_outfinal.fastq "$sample"_R2_trim_outfinal.fastq
    samtools sort -o /scratch/projects/sarg2022/Seawater-Sterivex/phage_abundances/$sample.aln.sort.bam /scratch/projects/sarg2022/Seawater-Sterivex/phage_abundances/$sample.aln.sam
    rm /scratch/projects/sarg2022/Seawater-Sterivex/phage_abundances/$sample.aln.sam
    samtools index /scratch/projects/sarg2022/Seawater-Sterivex/phage_abundances/$sample.aln.sort.bam
    samtools idxstats /scratch/projects/sarg2022/Seawater-Sterivex/phage_abundances/$sample.aln.sort.bam > /scratch/projects/sarg2022/Seawater-Sterivex/phage_abundances/$sample.idxstats
  done

  # Map the SSVdb to the Seawater Viral Fraction reads 
  #!/bin/bash
  eval "$(conda shell.bash hook)"
  conda activate /nethome/pxh390/miniconda3/envs/samtools
  for f in *_R1_001_outfinal.fastq; do
    sample=$(echo ${f} | sed 's/_R1_001_outfinal.fastq//g')
    smalt map -n 40 -y 0.95 -o /scratch/projects/sarg2022/Seawater-virome/phage_abundances/$sample.aln.sam \
    /scratch/projects/sarg2022/vMAG_dereplication/all_phages_050523 \
    "$sample".R1_001_outfinal.fastq "$sample".R2_001_outfinal.fastq
    samtools sort -o /scratch/projects/sarg2022/Seawater-virome/phage_abundances/$sample.aln.sort.bam /scratch/projects/sarg2022/Seawater-virome/phage_abundances/$sample.aln.sam
    rm /scratch/projects/sarg2022/Seawater-virome/phage_abundances/$sample.aln.sam
    samtools index /scratch/projects/sarg2022/Seawater-virome/phage_abundances/$sample.aln.sort.bam
    samtools idxstats /scratch/projects/sarg2022/Seawater-virome/phage_abundances/$sample.aln.sort.bam > /scratch/projects/sarg2022/Seawater-virome/phage_abundances/$sample.idxstats
  done
#------------------------------------------------- Creation of viral proteomic tree with GL-UVAB ---------------------------------------------------#
# Identify viruses of closest relatives in the ICTV database
#!/bin/bash 
eval "$(conda shell.bash hook)" 
conda activate gluvabkhf
perl /nethome/nsv19/anaconda3/envs/gluvab/GLUVAB_polyN_v0.6.pl \
--genomes_file_1 /projectnb/viralecology/databases/NCBI/viral_genomes_102422/single_line_95_perc_all_phage_genomes_of_bacteria_and_archaea_ICTV_102522.fasta \
--genomes_file_2 /scratch/projects/sarg2022/vMAG_dereplication/ALL_DEREPLICATED_vMAGs_great10000_050523.fasta
#
awk '{print $1}' GLUVAB_Scaffold_Recip_Scores.tsv | sort | uniq > phageIDs4ictv1_dereplicated.txt
awk '{print $2}' GLUVAB_Scaffold_Recip_Scores.tsv | sort | uniq > phageIDs4ictv2_dereplicated.txt

conda deactivate 

conda activate /nethome/pxh390/miniconda3/envs/utils
seqkit grep -n -f phageIDs4ictv1_dereplicated.txt ALL_DEREPLICATED_vMAGs_10000_030923.fasta > ictv_phages1.fna
seqkit grep -n -f phageIDs4ictv2_dereplicated.txt /projectnb/viralecology/databases/NCBI/viral_genomes_102422/single_line_95_perc_all_phage_genomes_of_bacteria_and_archaea_ICTV_102522.fasta > ictv_phages2.fna
cat ictv_phages1.fna ictv_phages2.fna > ictv_final_seqs_032823.fnac

conda deactivate 

# Identify viruses of closest relatives in the GLUVAB database
# do it with gluvab database
#!/bin/bash 
eval "$(conda shell.bash hook)" 
conda activate gluvabkhf
perl /nethome/nsv19/anaconda3/envs/gluvab/GLUVAB_polyN_v0.6.pl \
--genomes_file_1 /projectnb/viralecology/databases/gluvab/single_line_98perc_filtered_Gluvab_Viromes_11-16-2020.fasta \
--genomes_file_2 /scratch/projects/sarg2022/vMAG_dereplication/ALL_DEREPLICATED_vMAGs_great10000_050523.fasta
#
awk '{print $1}' GLUVAB_Scaffold_Recip_Scores.tsv | sort | uniq > phageIDs4gluvab1_cdhit.txt
awk '{print $2}' GLUVAB_Scaffold_Recip_Scores.tsv | sort | uniq > phageIDs4gluvab2_cdhit.txt

# For matches to the gluvab database I only used matches that had host information to the class level 
seqkit grep -n -f phageIDs4gluvab2_aquatic_marine_bacteria_class.txt \
/projectnb/viralecology/databases/gluvab/single_line_98perc_filtered_Gluvab_Viromes_11-16-2020.fasta \
> gluvab_phages2_aquatic_marine_bacteria_class.fna

conda deactivate 
# combine all the files into one file 

# Re-run gluvab to generate the tree 
#!/bin/bash 
eval "$(conda shell.bash hook)" 
conda activate gluvabkhf
perl /nethome/nsv19/anaconda3/envs/gluvab/GLUVAB_polyN_v0.6.pl --genomes_file_1 gluvab_ictv_class_final_seq.fna
conda deactivate

#-------------------- Identification of potential temperate viruses with a curated database of temperate phage related proteins  ------------------------#
# I manually curated a database of 2916 temperate phage related proteins consisting of integrases, recombinases, excionsaes, and transposases 
# from swissprot, pVOG, NCBI's viral protein database and ran a BLASTp against the SSVdb

#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate utils
makeblastdb -in temperate_proteins_061423.faa -dbtype prot -out temperate_phage_protDB
blastp -db temperate_phage_protDB \
        -query /scratch/projects/sarg2022/vMAG_dereplication/VIBRANT_ALL_DEREPLICATED_vMAGs_great10000_050523/VIBRANT_phages_ALL_DEREPLICATED_vMAGs_great10000_050523/ALL_DEREPLICATED_vMAGs_great10000_050523.ph$
        -out ALL_DEREPLICATED_vMAGs_temp_prot.out -evalue 1e-5 \
        -outfmt "6 qseqid sseqid pident qcovs qlen slen length mismatch qstart qend sstart send evalue" -num_threads 16
conda deactivate