#!/usr/bin/env python2

import sys, os, subprocess, multiprocessing, logging, shutil, gzip, re, ftplib
from StringIO import StringIO
from Bio import SeqIO
from collections import defaultdict

class Virseaqc():
	def __init__(self, threads, mem, db, reads_dir, outdir, verbosity, host_genome, adapters):
		self.mem = str(mem)
		self.db = db
		self.reads_dir = reads_dir
		self.reads_list = [ os.path.join(reads_dir, reads) for reads in os.listdir(reads_dir) if (os.path.isfile(os.path.join(reads_dir, reads)) and reads.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')))]
		self.max_threads = str(threads)
		self.application_threads = str(max(min(4, threads), threads/len(self.reads_list))) #Thanks T. Seeman
		self.parallel_job_limit = str(min(len(self.reads_list), threads/int(self.application_threads)))
		logging.info("Optimizing number of threads per parallel job...")
		logging.info("Using {} threads in {} jobs simultaneously.".format(self.application_threads, self.parallel_job_limit))
		self.outdir = outdir
		self.intermediate_files = os.path.join(outdir, 'intermediate_files')
		# self.fastqc_dir = os.path.join(outdir, 'fastqc_reads_quality_control')
		self.reference_files = os.path.join(self.intermediate_files, 'reference_sequence_files')
		self.databases = os.path.join(self.intermediate_files, 'databases')
		self.host_filtered_reads = os.path.join(self.intermediate_files, 'host_filtered_reads')
		# self.adapter_trimmed_dir = os.path.join(outdir, 'host_filtered_and_adapter_trimmed_reads')
		self.raw_ass_dir = os.path.join(self.intermediate_files, 'raw_assembly_output')
		# self.final_ass_dir = os.path.join(outdir, 'assemblies')
		self.blast_dir = os.path.join(self.intermediate_files, 'raw_blast_output')
		# self.blobtools_dir = os.path.join(outdir, 'blast_results')
		self.temp_dir = os.path.join(self.intermediate_files, '.temp')
		self.verbosity = verbosity
		self.host_genome = host_genome
		self.host_genome_index = os.path.join(self.reference_files, 'host_genome_index')
		self.adapters = adapters
	
	def clear_temp_dir(self):
		if os.path.isdir(self.temp_dir):
			shutil.rmtree(self.temp_dir)
		os.makedirs(self.temp_dir)
	
	def build_host_db(self):
		make_outdirs(self.reference_files)
		outputs = [
					'{}.1.bt2'.format(self.host_genome_index), 
					'{}.2.bt2'.format(self.host_genome_index), 
					'{}.3.bt2'.format(self.host_genome_index), 
					'{}.4.bt2'.format(self.host_genome_index), 
					'{}.rev.1.bt2'.format(self.host_genome_index), 
					'{}.rev.2.bt2'.format(self.host_genome_index)
					]
		if not check_multiple_outputs(outputs):
			logging.info("Building host genome bowtie2 index...")
			exitcode = run_command([
									'bowtie2-build',
									'--threads', self.max_threads,
									self.host_genome,
									self.host_genome_index
									])
			if exitcode == 0:
				logging.info("Host genome index built successfully.")
			else:
				logging.error("Execution of bowtie2-build failed. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
		else:
			logging.info("Host genome index detected. Skipping host genome index build and continuing to next step...")
	
	def build_viral_db(self):
		make_outdirs(self.databases)
		logging.info("Building viral genome blast database...")
		download_multipart_db('refseq/release/viral/', 'viral.{}.protein.faa.gz')
		with open(os.path.join(self.databases, 'viral.protein.faa'), 'wb') as outfile:
			for file in os.listdir('.'):
				if file.endswith('.protein.faa.gz') and file.startswith('viral.'):
					with gzip.open(file, 'rb') as infile:
						shutil.copyfileobj(infile, outfile)
					os.remove(file)
		exitcode = run_command(['diamond', 'makedb', '--in', os.path.join(self.databases, 'viral.protein.faa'), '-d', os.path.join(self.databases, 'viraldb')])
		if exitcode == 0:
			logging.info("Viral genome blast database built successfully.")
		else:
			logging.error("Viral genome blast database building failed. See logfile virseaqc.log in your output directory for more details.")
			sys.exit(1)
	
	def filter_host_reads(self):
		make_outdirs(self.host_filtered_reads)
		for reads in self.reads_list:
			make_outdirs(os.path.join(self.outdir, get_sample_name(reads)))
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.host_filtered_reads, "{}.host_filtered.fastq.gz".format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without host read filtering. Starting bowtie2...".format(len(not_done)))
			for reads in not_done:
				logging.info("Filtering host reads from sample {}...".format(get_sample_name(reads)))
				exitcode = run_command([
										'virseaqc_filter_host_reads',
										reads,
										self.host_genome_index,
										os.path.join(self.host_filtered_reads, get_sample_name(reads)),
										self.max_threads
										])
				if exitcode == 0:
					logging.info("Host read filtering completed successfully.")
				else:
					logging.error("Host read filtering failed. See logfile virseaqc.log in your output directory for more details.")
					sys.exit(1)
		elif len(not_done) == 0 and len(self.reads_list) > 0:
			logging.info("Host filtered reads detected for all samples. Skipping host read filtering and continuing to next step...")
		else:
			logging.info("No reads found in {}. Please check that the -q|--fastq-dir option is correct.".format(self.reads_dir))
			sys.exit(1)
	
	def trim_adapters(self):
		# make_outdirs(self.adapter_trimmed_dir)
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.outdir, get_sample_name(reads), "{}.host_filtered_adapter_trimmed.fastq.gz".format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without adapter trimming. Starting cutadapt...".format(len(not_done)))
			cmds = [ ' '.join([
					'cutadapt',
					'-a', 'file:{}'.format(self.adapters),
					'-o', os.path.join(self.outdir, get_sample_name(reads), "{}.host_filtered_adapter_trimmed.fastq.gz".format(get_sample_name(reads))),
					os.path.join(self.host_filtered_reads, "{}.host_filtered.fastq.gz".format(get_sample_name(reads)))
					]) 
					for reads in not_done 
					]
			exitcode = run_command(['parallel', '-j{}'.format(self.max_threads), ':::'] + cmds)
			if exitcode == 0:
				logging.info("Adapter removal completed successfully.")
			else:
				logging.error("Adapter removal failed. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
		elif len(not_done) == 0 and len(self.reads_list) > 0:
			logging.info("Adapter removed reads detected for all samples. Skipping Adapter removal and continuing to next step...")
		else:
			logging.info("No reads found in {}. Please check that the -q|--fastq-dir option is correct.".format(self.reads_dir))
			sys.exit(1)
	
	def run_fastqc(self):
		# make_outdirs(self.fastqc_dir)
		trimmed_not_done = [ os.path.join(self.outdir, get_sample_name(reads), "{}.host_filtered_adapter_trimmed.fastq.gz".format(get_sample_name(reads))) for reads in self.reads_list if not os.path.isfile(os.path.join(self.outdir, get_sample_name(reads), "{}.host_filtered_adapter_trimmed_fastqc.html".format(get_sample_name(reads)))) ]
		if len(trimmed_not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without fastqc reports. Starting fastqc...".format(len(trimmed_not_done))) 
			cmds = [ ' '.join([
								'fastqc',
								'--outdir={}'.format(os.path.join(self.outdir, get_sample_name(reads))),
								reads
								]) 
					for reads in trimmed_not_done 
					]
			exitcode = run_command(['parallel', '-j{}'.format(self.max_threads), ':::'] + cmds)
			if exitcode == 0:
				logging.info("fastqc reports generated successfully.")
			else:
				logging.error("fastqc failed. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
		elif len(trimmed_not_done) == 0 and len(self.reads_list) > 0:
			logging.info("fastqc reports detected for all samples. Skipping fastqc and continuing to next step...")
		else:
			logging.info("No reads found in {}. Please check that the -q|--fastq-dir option is correct.".format(self.reads_dir))
			sys.exit(1)
	
	def run_diamond(self):
		make_outdirs(self.blast_dir)
		trimmed_not_done = [ os.path.join(self.outdir, get_sample_name(reads), "{}.host_filtered_adapter_trimmed.fastq.gz".format(get_sample_name(reads))) for reads in self.reads_list if not os.path.isfile(os.path.join(self.outdir, get_sample_name(reads), '{}.virus_classification.html'.format(get_sample_name(reads)))) ]
		print(trimmed_not_done)
		if len(trimmed_not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without viral ID. Starting diamond...".format(len(trimmed_not_done)))
			for reads in trimmed_not_done:
				logging.info("Identifying viral reads from sample {}...".format(get_sample_name(reads)))
				fasta_path = convert_fastq_to_fasta(reads, self.blast_dir)
				outfile_path = os.path.join(self.blast_dir, '{}.daa'.format(get_sample_name(reads)))
				outfile_blast = os.path.join(self.blast_dir, '{}.dmnd.blast'.format(get_sample_name(reads)))
				outfile_krona = os.path.join(self.outdir, get_sample_name(reads), '{}.virus_classification.html'.format(get_sample_name(reads)))
				exitcode = run_command([
										'diamond',
										'blastx',
										'-d', os.path.join(self.databases, 'viraldb'),
										'-f', '100',
										'--threads', self.max_threads,
										'--salltitles',
										'-k', '1',
										'-q', fasta_path,
										'-o', outfile_path
										])
				if exitcode != 0:
					logging.error("Diamond BLAST failed. See logfile virseaqc.log in your output directory for more details.")
					sys.exit(1)
				exitcode = run_command([
										'diamond',
										'view',
										'-a', outfile_path,
										'-o', outfile_blast
										])
				if exitcode != 0:
					logging.error("Diamond BLAST file conversion failed. See logfile virseaqc.log in your output directory for more details.")
					sys.exit(1)
				exitcode = run_command([
										'ktImportBLAST',
										'-i',
										'-p',
										'-o', outfile_krona,
										outfile_blast
										])
				if exitcode == 0:
					logging.info("Viral ID successful.")
				else:
					logging.error("Viral ID failed. See logfile virseaqc.log in your output directory for more details.")
					sys.exit(1)
		elif len(trimmed_not_done) == 0 and len(self.reads_list) > 0:
			logging.info("Viral IDs detected for all samples. Skipping diamond and continuing to next step...")
		else:
			logging.info("No reads found in {}. Please check that the -q|--fastq-dir option is correct.".format(self.reads_dir))
			sys.exit(1)
	

	def assemble_reads(self):
		make_outdirs(self.raw_ass_dir)
		trimmed_not_done = [ os.path.join(self.outdir, get_sample_name(reads), "{}.host_filtered_adapter_trimmed.fastq.gz".format(get_sample_name(reads))) for reads in self.reads_list if not os.path.isfile(os.path.join(self.raw_ass_dir, get_sample_name(reads), 'scaffolds.fasta')) ]
		if len(trimmed_not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without assemblies. Starting SPAdes...".format(len(trimmed_not_done)))
			null = ( shutil.rmtree(os.path.join(self.raw_ass_dir, get_sample_name(reads))) for reads in trimmed_not_done if os.path.isdir(os.path.join(self.raw_ass_dir, get_sample_name(reads))) )
			cmds = [ ' '.join([
						'spades.py',
						'--careful',
						'--iontorrent',
						'-s', reads,
						'-t', self.application_threads,
						'-m', self.mem,
						'-k', kmer_calc(reads),
						'--cov-cutoff', 'auto',
						'-o', os.path.join(self.raw_ass_dir, get_sample_name(reads))
						]) 
						for reads in trimmed_not_done ]
			exitcode = run_command(['parallel', '-j{}'.format(self.parallel_job_limit), ':::'] + cmds)
			trimmed_not_done_again = [ os.path.join(self.outdir, get_sample_name(reads), "{}.host_filtered_adapter_trimmed.fastq.gz".format(get_sample_name(reads))) for reads in self.reads_list if not os.path.isfile(os.path.join(self.raw_ass_dir, get_sample_name(reads), 'scaffolds.fasta')) ]
			if len(trimmed_not_done_again) > 0:
				logging.info("{} reads did not assemble using SPAdes with the --careful option. Retrying SPAdes assemblies without --careful option...".format(len(trimmed_not_done_again)))
				null = ( shutil.rmtree(os.path.join(self.raw_ass_dir, get_sample_name(reads))) for reads in trimmed_not_done_again if os.path.isdir(os.path.join(self.raw_ass_dir, get_sample_name(reads))) )
				cmds = [ ' '.join([
							'spades.py',
							'--iontorrent',
							'-s', reads,
							'-t', self.application_threads,
							'-m', self.mem,
							'-k', kmer_calc(reads),
							'--cov-cutoff', 'auto',
							'-o', os.path.join(self.raw_ass_dir, get_sample_name(reads))
							] )
							for reads in trimmed_not_done_again ]
				exitcode = run_command(['parallel', '-j{}'.format(self.parallel_job_limit), ':::'] + cmds)
			if exitcode == 0:
				logging.info("Genome assembly completed successfully for all samples.")
				# make_outdirs(self.final_ass_dir)
				for reads in self.reads_list:
					shutil.copy2(os.path.join(self.raw_ass_dir, get_sample_name(reads), 'scaffolds.fasta'), os.path.join(self.outdir, get_sample_name(reads), '{}.scaffolds.fasta'.format(get_sample_name(reads))))
			else:
				logging.error("Genome assembly failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
		elif len(trimmed_not_done) == 0 and len(self.reads_list) > 0:
			logging.info("Assemblies detected for all samples. Skipping spades assembly and continuing to next step...")
		else:
			logging.info("No reads found in {}. Please check that the -q|--fastq-dir option is correct.".format(self.reads_dir))
			sys.exit(1)
		
	def run_blastn(self):
		make_outdirs(self.blast_dir)
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.blast_dir, '{}.blast_hits'.format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Searching nt database with blastn...")
			for reads in not_done:
				exitcode = run_command([
										'blastn', 
										'-query', os.path.join(self.raw_ass_dir, get_sample_name(reads), 'scaffolds.fasta'),
										'-db', self.db,
										'-outfmt', '6 qseqid staxids bitscore std',
										'-max_target_seqs', '10',
										'-max_hsps', '1',
										'-evalue', '1e-25',
										'-num_threads', self.max_threads,
										'-out', os.path.join(self.blast_dir, '{}.blast_hits'.format(get_sample_name(reads)))
										])
				if exitcode == 0:
					logging.info("BLAST searches completed successfully for sample {}.".format(get_sample_name(reads)))
				else:
					logging.error("BLAST searches failed for sample {}. See logfile virseaqc.log in your output directory for more details.".format(get_sample_name(reads)))
					sys.exit(1)
	
	def run_bt_create(self):
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.blast_dir, '{}.blobDB.json'.format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Running blobtools for summaries and graphical output...")
			cmds = [ ' '.join([
						'blobtools', 'create',
						'-i', os.path.join(self.raw_ass_dir, get_sample_name(reads), 'scaffolds.fasta'), 
						'-y', 'spades',
						'-t', os.path.join(self.blast_dir, '{}.blast_hits'.format(get_sample_name(reads))),
						'-x', 'bestsumorder',
						'-o', os.path.join(self.blast_dir, get_sample_name(reads))
						]) for reads in self.reads_list ]
			exitcode = run_command(['parallel', '-j{}'.format(self.max_threads), ':::'] + cmds)
			if exitcode != 0:
				logging.error("blobtools failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
	
	def run_bt_view(self):
		# make_outdirs(self.blobtools_dir)
		# for reads in self.reads_list:
			# make_outdirs(os.path.join(self.blobtools_dir, get_sample_name(reads)))
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.outdir, get_sample_name(reads), '{}.blast_results.tsv'.format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			cmds = [ ' '.join([
						'blobtools', 'view',
						'-i', os.path.join(self.blast_dir, '{}.blobDB.json'.format(get_sample_name(reads))),
						'-x', 'bestsumorder',
						'-r', 'species',
						'-b',
						'-o', os.path.join(self.blast_dir, get_sample_name(reads))
						]) for reads in self.reads_list ]
			exitcode = run_command(['parallel', '-j{}'.format(self.max_threads), ':::'] + cmds)
			if exitcode != 0:
				logging.error("blobtools failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
			else:
				for reads in self.reads_list:
					shutil.copy2(os.path.join(self.blast_dir, '{0}.{0}.blobDB.table.txt'.format(get_sample_name(reads))), os.path.join(self.outdir, get_sample_name(reads), '{}.blast_results.tsv'.format(get_sample_name(reads))))
	
	def run_bt_plot_genus(self):
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.outdir, get_sample_name(reads), '{}.blast_results.genus_blobplot.png'.format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			cmds = [ ' '.join([
						'blobtools', 'plot',
						'-i', os.path.join(self.blast_dir, '{}.blobDB.json'.format(get_sample_name(reads))),
						'-x', 'bestsumorder',
						'-r', 'genus',
						'-l', '100',
						'-o', os.path.join(self.blast_dir, get_sample_name(reads))
						]) for reads in self.reads_list ]
			exitcode = run_command(['parallel', '-j{}'.format(self.max_threads), ':::'] + cmds)
			if exitcode != 0:
				logging.error("blobtools failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
			else:
				for reads in self.reads_list:
					shutil.copy2(os.path.join(self.blast_dir, '{0}.{0}.blobDB.json.bestsumorder.genus.p7.span.100.blobplot.spades.png'.format(get_sample_name(reads))), os.path.join(self.outdir, get_sample_name(reads), '{}.blast_results.genus_blobplot.png'.format(get_sample_name(reads))))

	
	def run_bt_plot_species(self):
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.outdir, get_sample_name(reads), '{}.blast_results.species_blobplot.png'.format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			cmds = [ ' '.join([
						'blobtools', 'plot',
						'-i', os.path.join(self.blast_dir, '{}.blobDB.json'.format(get_sample_name(reads))),
						'-x', 'bestsumorder',
						'-r', 'species',
						'-l', '100',
						'-o', os.path.join(self.blast_dir, get_sample_name(reads))
						]) for reads in self.reads_list ]
			exitcode = run_command(['parallel', '-j{}'.format(self.max_threads), ':::'] + cmds)
			if exitcode != 0:
				logging.error("blobtools failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
			else:
				for reads in self.reads_list:
					shutil.copy2(os.path.join(self.blast_dir, '{0}.{0}.blobDB.json.bestsumorder.species.p7.span.100.blobplot.spades.png'.format(get_sample_name(reads))), os.path.join(self.outdir, get_sample_name(reads), '{}.blast_results.species_blobplot.png'.format(get_sample_name(reads))))

	
	def run_pipeline(self):
		self.build_host_db()
		self.build_viral_db()
		self.filter_host_reads()
		self.trim_adapters()
		self.run_fastqc()
		self.run_diamond()
		self.assemble_reads()
		self.run_blastn()
		self.run_bt_create()
		self.run_bt_view()
		self.run_bt_plot_genus()
		self.run_bt_plot_species()



def open_reads(reads):
	if reads.endswith('.gz'):
		reads_handle = gzip.open(reads, 'rt')
	else:
		reads_handle = open(reads, 'r')
	return reads_handle

def kmer_calc(reads):
	sample_name = get_sample_name(reads)
	reads_handle = open_reads(reads)
	max_read_len = max(len(line) for i, line in enumerate(reads_handle) if (i % 4) == 1 )
	reads_handle.close()
	run_ks = []
	for k in [21,33,55,77,99,127]:
		if k < (max_read_len/2):
			run_ks.append(str(k))
	logging.info("{}: Running spades with kmers {}".format(sample_name, ",".join(run_ks)))
	return ','.join(run_ks)

def check_multiple_outputs(outputs):
	for output in outputs:
		if not os.path.isfile(output):
			return False
	return True

def get_sample_name(reads):
	name, suffix = os.path.splitext(os.path.basename(reads))
	if suffix == '.gz':
		name, suffix = os.path.splitext(name)
	name2, suffix = os.path.splitext(os.path.basename(name))
	if suffix == '.host_filtered' or suffix == '.host_filtered_adapter_trimmed':
		return name2
	return name


def log_subprocess_output(pipe):
	for line in iter(pipe.readline, b''): # b'\n'-separated lines
		logging.debug(line.strip())


def make_outdirs(dir):
	if os.path.isdir(dir) == False:
		os.makedirs(dir)


def run_command(cmd, shell_opt=False):
	if shell_opt == True:
		logging.info("Executing command: {}".format(cmd))
	else:
		logging.info("Executing command: {}".format(" ".join(cmd)))
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=shell_opt)
	with process.stdout:
		log_subprocess_output(process.stdout)
	exitcode = process.wait()
	return exitcode


def configure_logs(verbosity, outdir):
	logging.basicConfig(filename=os.path.join(outdir, 'virseaqc.log'),format='%(levelname)s:%(message)s',level=logging.DEBUG)
	cmd_log = logging.StreamHandler()
	if verbosity == True:
		cmd_log.setLevel(logging.DEBUG)
	else:
		cmd_log.setLevel(logging.INFO)
	logging.getLogger().addHandler(cmd_log)


def configure_outdir(outdir, force):
	if os.path.isdir(outdir):
		if force == True:
			shutil.rmtree(outdir)
			os.makedirs(outdir)
	else:
		os.makedirs(outdir)

def download_multipart_db(path, filename):
	ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov')
	ftp.login("anonymous", "daniel.bogema@dpi.nsw.gov.au")
	ftp.cwd(path)
	file_num = 0
	while True:
		try:
			file_num += 1
			file = filename.format(file_num)
			ftp.retrbinary("RETR " + file, open(file, 'wb').write)
		except ftplib.error_perm:
			os.remove(file)
			break

def convert_fastq_to_fasta(fastq_file, output_dir):
	sample_name = get_sample_name(fastq_file)
	recs = SeqIO.parse(open_reads(fastq_file), 'fastq')
	fasta_path = os.path.join(output_dir, '{}.fasta'.format(sample_name))
	null = SeqIO.write(recs, fasta_path, 'fasta')
	return fasta_path


def run(threads, mem, reads_dir, force, outdir, verbosity, db, host_genome, adapters):
	configure_outdir(outdir, force)
	configure_logs(verbosity, outdir)
	logging.info("Beginning virseaqc run...")
	job = Virseaqc(threads, mem, db, reads_dir, outdir, verbosity, host_genome, adapters)
	if force == True:
		logging.info("Force option given. Deleting previous run...")
	logging.info("Log for the current run can be found at {}.".format(os.path.join(outdir, 'virseaqc.log')))
	job.run_pipeline()
	logging.info("Run complete. Thanks for using virseaqc.")
	job.clear_temp_dir()

