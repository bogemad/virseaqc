#!/usr/bin/env python2

import sys, os, subprocess, multiprocessing, logging, shutil, gzip, re
from StringIO import StringIO
from collections import defaultdict

class Virseaqc():
	def __init__(self, threads, mem, db, reads_dir, outdir, verbosity, host_genome, adapters):
		self.mem = str(mem)
		self.db = db
		self.reads_dir = reads_dir
		self.reads_list = [ os.path.join(reads_dir, reads) for reads in os.listdir(reads_dir) if (os.path.isfile(os.path.join(reads_dir, reads)) and reads.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')))]
		self.threads = str(threads if threads < len(self.reads_list) else len(self.reads_list))
		self.outdir = outdir
		self.fastqc_dir = os.path.join(outdir, 'qc_plots', 'fastqc')
		self.reference_files = os.path.join(outdir, 'reference_sequence_files')
		self.host_filtered_reads = os.path.join(outdir, 'host_filtered_reads')
		self.adapter_trimmed_dir = os.path.join(outdir, 'host_filtered_and_adapter_trimmed_reads')
		self.raw_ass_dir = os.path.join(outdir, 'raw_assemblies')
		self.blobtools_dir = os.path.join(outdir, 'qc_plots', 'blobtools')
		self.temp_dir = os.path.join(outdir, '.temp')
		self.successful = defaultdict(list)
		self.failed = defaultdict(list)
		self.verbosity = verbosity
		self.host_genome = host_genome
		self.host_genome_index = os.path.join(self.reference_files, 'host_genome_index')
		self.adapters = adapters
	
	def clear_temp_dir(self):
		if os.path.isdir(self.temp_dir):
			shutil.rmtree(self.temp_dir)
		os.makedirs(self.temp_dir)
	
	def build_host_db(self):
		sample_name = get_sample_name(reads)
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
									'--threads', self.threads,
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
	
	def filter_host_reads(self):
		make_outdirs(self.host_filtered_reads)
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.host_filtered_reads, "{}.host_filtered.fastq".format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without host read filtering. Starting bowtie2...".format(len(not_done)))
			cmds = [ [
					'virseaqc_filter_host_reads',
					reads,
					self.host_genome_index,
					get_sample_name(reads)
					] for reads in not_done 
					]
			exitcode = run_command(['parallel', '-j{}'.format(self.threads), ':::'] + cmds)
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
		make_outdirs(self.adapter_trimmed_dir)
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.adapter_trimmed_dir, "{}.host_filtered_adapter_trimmed.fastq".format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without adapter trimming. Starting cutadapt...".format(len(not_done)))
			cmds = [ [
					'cutadapt',
					'-a', 'file:{}'.format(self.adapters),
					'-o', os.path.join(self.adapter_trimmed_dir, "{}.host_filtered_adapter_trimmed.fastq".format(get_sample_name(reads)))
					os.path.join(self.host_filtered_reads, "{}.host_filtered.fastq".format(get_sample_name(reads)))
					] 
					for reads in not_done 
					]
			exitcode = run_command(['parallel', '-j{}'.format(self.threads), ':::'] + cmds)
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
		make_outdirs(self.fastqc_dir)
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.fastqc_dir, "{}_fastqc.html".format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without fastqc reports. Starting fastqc...".format(len(not_done))) 
			exitcode = run_command([
									'parallel', '-j{}'.format(self.threads),
									'fastqc',
									'--outdir={}'.format(self.fastqc_dir),
									'{}', ':::'] + not_done)
			if exitcode == 0:
				logging.info("fastqc reports generated successfully.")
			else:
				logging.error("fastqc failed. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
		elif len(not_done) == 0 and len(self.reads_list) > 0:
			logging.info("fastqc reports detected for all samples. Skipping fastqc and continuing to next step...")
		else:
			logging.info("No reads found in {}. Please check that the -q|--fastq-dir option is correct.".format(self.reads_dir))
			sys.exit(1)
	
	def assemble_reads(self):
		make_outdirs(self.raw_ass_dir)
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.raw_ass_dir, get_sample_name(reads), 'scaffolds.fasta')) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Found {} reads without assemblies. Starting SPAdes...".format(len(not_done)))
			null = ( shutil.rmtree(os.path.join(self.raw_ass_dir, get_sample_name(reads))) for reads in not_done if os.path.isdir(os.path.join(self.raw_ass_dir, get_sample_name(reads))) )
			cmds = [ [
						'spades.py',
						'--careful',
						'--se', reads,
						'-t', '1',
						'-m', self.mem,
						'-k', kmer_calc(reads),
						'--cov-cutoff', 'auto',
						'-o', os.path.join(self.raw_ass_dir, get_sample_name(reads))
						] 
						for reads in not_done ]
			exitcode = run_command(['parallel', '-j{}'.format(self.threads), ':::'] + cmds)
			not_done_again = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.raw_ass_dir, get_sample_name(reads), 'scaffolds.fasta')) ]
			if not_done_again > 0:
				logging.info("{} reads did not assemble using SPAdes with the --careful option. Retrying SPAdes assemblies without --careful option...".format(len(not_done_again)))
				null = ( shutil.rmtree(os.path.join(self.raw_ass_dir, get_sample_name(reads))) for reads in not_done if os.path.isdir(os.path.join(self.raw_ass_dir, get_sample_name(reads))) )
				cmds = [ [
							'spades.py',
							'--se', reads,
							'-t', '1',
							'-m', self.mem,
							'-k', kmer_calc(reads),
							'--cov-cutoff', 'auto',
							'-o', os.path.join(self.raw_ass_dir, get_sample_name(reads))
							] 
							for reads in not_done ]
				exitcode = run_command(['parallel', '-j{}'.format(self.threads), ':::'] + cmds)
			if exitcode == 0:
				logging.info("Genome assembly completed successfully for all samples.")
			else:
				logging.error("Genome assembly failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
		
	def run_blastn(self):
		null = (make_outdirs(os.path.join(self.blobtools_dir, get_sample_name(reads))) for reads in self.reads_list)
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.blast_hits'.format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Searching nt database with blastn...")
			cmds = [ [
						'blastn', 
						'-query', os.path.join(self.raw_ass_dir, get_sample_name(reads), 'scaffolds.fasta'),
						'-db', self.db,
						'-outfmt', '6 qseqid staxids bitscore std',
						'-max_target_seqs', '10',
						'-max_hsps', '1',
						'-evalue', '1e-25',
						'-num_threads', '1',
						'-out', os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.blast_hits'.format(get_sample_name(reads)))
						] for reads in self.reads_list ]
			exitcode = run_command(['parallel', '-j{}'.format(self.threads), ':::'] + cmds)
			if exitcode == 0:
				logging.info("BLAST searches completed successfully for all samples.")
			else:
				logging.error("BLAST searches failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
	
	def run_bt_create(self, reads):
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.blobDB.json'.format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			logging.info("Running blobtools for summaries and graphical output...")
			cmds = [ [
						'blobtools', 'create',
						'-i', os.path.join(self.raw_ass_dir, get_sample_name(reads), 'scaffolds.fasta'), 
						'-y', 'spades',
						'-t', os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.blast_hits'.format(get_sample_name(reads))),
						'-x', 'bestsumorder',
						'-o', os.path.join(self.blobtools_dir, get_sample_name(reads), get_sample_name(reads))
						] for reads in self.reads_list ]
			exitcode = run_command(['parallel', '-j{}'.format(self.threads), ':::'] + cmds)
			if exitcode != 0:
				logging.error("blobtools failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
	
	def run_bt_view(self, reads):
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.blobDB.bestsumorder.table.txt'.format(get_sample_name(reads)))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			cmds = [ [
						'blobtools', 'view',
						'-i', os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.blobDB.json'.format(get_sample_name(reads))),
						'-x', 'bestsumorder',
						'-r', 'species',
						'-b',
						'-o', os.path.join(self.blobtools_dir, get_sample_name(reads), get_sample_name(reads))
						] for reads in self.reads_list ]
			exitcode = run_command(['parallel', '-j{}'.format(self.threads), ':::'] + cmds)
			if exitcode != 0:
				logging.error("blobtools failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
	
	def run_bt_plot_genus(self, reads):
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.bestsumorder.genus.p7.span.100.blobplot.spades.png'.format(blobdb_file))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			cmds = [ [
						'blobtools', 'plot',
						'-i', os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.blobDB.json'.format(get_sample_name(reads))),
						'-x', 'bestsumorder',
						'-r', 'genus',
						'-l', '100',
						'-o', os.path.join(self.blobtools_dir, get_sample_name(reads), get_sample_name(reads))
						] for reads in self.reads_list ]
			exitcode = run_command(['parallel', '-j{}'.format(self.threads), ':::'] + cmds)
			if exitcode != 0:
				logging.error("blobtools failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
	
	def run_bt_plot_species(self, reads):
		not_done = [ reads for reads in self.reads_list if not os.path.isfile(os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.bestsumorder.species.p7.span.100.blobplot.spades.png'.format(blobdb_file))) ]
		if len(not_done) > 0 and len(self.reads_list) > 0:
			cmds = [ [
						'blobtools', 'plot',
						'-i', os.path.join(self.blobtools_dir, get_sample_name(reads), '{}.blobDB.json'.format(get_sample_name(reads))),
						'-x', 'bestsumorder',
						'-r', 'genus',
						'-l', '100',
						'-o', os.path.join(self.blobtools_dir, get_sample_name(reads), get_sample_name(reads))
						] for reads in self.reads_list ]
			exitcode = run_command(['parallel', '-j{}'.format(self.threads), ':::'] + cmds)
			if exitcode != 0:
				logging.error("blobtools failed for one or more samples. See logfile virseaqc.log in your output directory for more details.")
				sys.exit(1)
	
	def run_pipeline(self):
		self.build_host_db()
		self.filter_host_reads()
		self.trim_adapters()
		self.run_fastqc()
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
		run_ks.append(str(k))
		if k > (max_read_len/2):
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
	return name


def log_subprocess_output(pipe):
	for line in iter(pipe.readline, b''): # b'\n'-separated lines
		logging.debug(line.strip())


def make_outdirs(dir):
	if os.path.isdir(dir) == False:
		os.makedirs(dir)


def run_command(cmd, shell_opt=False):
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

def run(threads, mem, reads_dir, force, outdir, verbosity, db):
	configure_outdir(outdir, force)
	configure_logs(verbosity, outdir)
	logging.info("Beginning virseaqc run...")
	job = Virseaqc(threads, mem, db, reads_dir, outdir, verbosity)
	if force == True:
		logging.info("Force option given. Deleting previous run...")
	logging.info("Log for the current run can be found at {}.".format(os.path.join(outdir, 'logs')))
	job.run_pipeline()
	logging.info("Run complete. Thanks for using virseaqc.")
	job.clear_temp_dir()

