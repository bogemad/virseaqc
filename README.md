# virseaqc
virseaqc - A combined pipeline for VIRal genome SEarching, Assembly, and Quality Control on Ion-Torrent instruments.

usage: virseaqc [-h] -q READS_DIR -o OUTDIR [-b DB] [-s HOST_GENOME]
                [-a ADAPTERS] [-t THREADS] [-p] [-i] [-m MEM] [-f] [--version]
                [-v]

virseaqc - A combined pipeline for VIRal genome SEarching, Assembly, and
Quality Control on Ion-Torrent instruments.

required arguments:
  -q READS_DIR, --fastq-dir READS_DIR
                        Path to a directory with your ion-torrent or single
                        end fastq files.
  -o OUTDIR, --output-directory OUTDIR
                        Path to the output directory. A directory will be
                        created if one does not exist.
  -b DB, --blast_database DB
                        Path to the local nt blast database. This pipeline
                        requires you to have a local copy of the nt database.
                        See https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&P
                        AGE_TYPE=BlastDocs&DOC_TYPE=Download to install a
                        local nt database.
  -s HOST_GENOME, --host_genome HOST_GENOME
                        Path to the host reference genome to remove host reads
  -a ADAPTERS, --adapters ADAPTERS
                        Path to the file containing adapter reads for removal
                        from sequence reads

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing.
  -p, --paired_end      Reads are sourced from an Illumina instrument rather
                        than an Ion Torrent and are paired-end. Note that
                        these files must end with n.fastq(.gz), where n=1,2
  -i, --interleaved     Reads are sourced from an Illumina instrument rather
                        than an Ion Torrent and are paired-end interleaved.
  -m MEM, --max_memory MEM
                        Maximum amount of RAM to assign to the pipeline in GB
                        (Just the number).
  -f, --force           Overwrite files in the output directories.
  --version             show program's version number and exit
  -v, --verbose         Increase verbosity on command line output (n.b.
                        verbose output is always saved to virseaqc.log in the
                        output directory).
