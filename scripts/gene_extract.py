#!/usr/bin/env python3
"""Extract gene specific data for a Nextstrain-TrepoGen gene dataset.

Description:
	This script prepares data for a TrepoGen gene dataset by
	i. extracting gene-specific sequences using MUSIAL,
	ii. processing the reference genome annotation to extract gene sub-features.

	The workflow first creates a MUSIAL storage file (including its configuration)
	and exports a non-reference sequence alignment for the target gene. It then
	parses the FASTA to extract and reformat the gene's reference sequence.
	Next, it reads the original GFF3, updates entries and coordinates — accounting
	for insertion-induced gaps — and writes a revised annotation plus a sub-feature
	list for augur translate. Finally, it translates specified sub-features into
	amino acid sequences, produces a CSV of per-sample sub-feature-type metadata,
	and generates an Auspice configuration for Nextstrain visualization.

Usage:
	python prepare_gene_build.py
		-ir <reference.fasta>
		-ia <annotation.gff3>
		-ig <gene_name>
		-iv <variants.vcf>
		-or <output_reference.fasta>
		-oa <output_annotation.gff3>
		-os <output_sequences.fasta>
		-og <output_genes.txt>
		-ot <output_types.csv>
		-oc <output_configuration.json>
		-m <musial.jar>
		[-q <root_sequence_identifier>]
"""
import json
import subprocess
import argparse
import os
import tempfile
from pandas import DataFrame
from numpy import linspace
from textwrap import wrap
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import colormaps, colors

"""Fixed biological region name map for gene sub-features."""
biological_region_name_map = { 
	'sp': 'Signal-peptide',
	'hatch': 'Hatch',
	'tm': 'Transmembrane-helix-',
	'bs': 'Beta-strand-',
	'icl': 'Intracellular-loop-',
	'ecl': 'Extracellular-loop-'
}

def parse_args():
	"""Parses command-line arguments for preparing data for a TrepoGen gene build using MUSIAL.

	Description:
		This function sets up the command-line interface for the script, allowing users to specify
		input files (reference genome, annotation, gene name, variants), output paths for processed
		data (reference sequence, annotation, sequences, gene names, region types metadata, and
		auspice configuration), and the path to the MUSIAL JAR file. It also allows an optional root
		sequence identifier to be specified.

	Arguments:
		-ir, --reference (str, required): Reference genome sequence (FASTA) file.
		-ia, --annotation (str, required): Reference genome annotation (GFF3) file.
		-ig, --gene (str, required): The gene to extract sequence and annotation files for.
		-iv, --vcf (str, required): Variants (VCF) file.
		-or, --output-reference (str, required): Output path for the reference sequence file (FASTA).
		-oa, --output-annotation (str, required): Output path for the reference annotation file (GFF3).
		-os, --output-sequences (str, required): Output path for the sequences file (FASTA).
		-og, --output-genes (str, required): Output path for the gene names file (TXT).
		-ot, --output-types (str, required): Output path for the per sample region type metadata (CSV).
		-oc, --output-configuration (str, required): Output path for the additional auspice configuration (JSON).
		-m, --musial (str, required): Path to the Musial JAR file or executable in PATH.
		-q, --root (str, optional): Optional sequence identifier to use for the root sequence.

	Returns:
		argparse.Namespace: Parsed command-line arguments with the following attributes:
		- reference (str): Path to the reference genome sequence file.
		- annotation (str): Path to the reference genome annotation file.
		- gene (str): Name of the gene to extract.
		- vcf (str): Path to the variants file.
		- output_reference (str): Path to the output reference sequence file.
		- output_annotation (str): Path to the output annotation file.
		- output_sequences (str): Path to the output sequences file.
		- output_genes (str): Path to the output gene names file.
		- output_types (str): Path to the output per sample region type metadata file.
		- output_configuration (str): Path to the output auspice configuration file.
		- musial (str): Path to the Musial JAR file.
		- root (str, optional): Optional sequence identifier for the root sequence.
	"""
	parser = argparse.ArgumentParser(
		description="Prepare data for a TrepoGen gene build using MUSIAL.",
	)
	parser.add_argument("-ir", "--reference", type=str, required=True, help="Reference genome sequence (fasta) file.")
	parser.add_argument("-ia", "--annotation", type=str, required=True, help="Reference genome annotation (gff3) file.")
	parser.add_argument("-ig", "--gene", type=str, required=True, help="The gene to extract sequence and annotation files for.")
	parser.add_argument("-iv", "--vcf", type=str, required=True, help="Variants (vcf) file.")
	parser.add_argument("-or", "--output-reference", type=str, required=True, help="Output path for the reference sequence file (fasta).")
	parser.add_argument("-oa", "--output-annotation", type=str, required=True, help="Output path for the reference annotation file (gff3).")
	parser.add_argument("-os", "--output-sequences", type=str, required=True, help="Output path for the sequences file (fasta).")
	parser.add_argument("-og", "--output-genes", type=str, required=True, help="Output path for the gene names file (txt).")
	parser.add_argument("-ot", "--output-types", type=str, required=True, help="Output path for the per sample ECL type metadata (csv).")
	parser.add_argument("-oc", "--output-configuration", type=str, required=True, help="Output path for the auspice configuration (json).")
	parser.add_argument("-m", "--musial", type=str, required=True, help="Path to the Musial jar file.")
	parser.add_argument("-q", "--root", type=str, required=False, help="Optional sequence identifier to use for the root sequence.")
	return parser.parse_args()

def write_config( reference : str, annotation : str, gene : str, vcf_input : str, tmpdir : str, min_coverage : int = 3, min_frequency : float = 0.8  ):
	"""Generates a MUSIAL configuration file for a Nextstrain-TrepoGen gene workflow.

	Description:
		This function creates a configuration file for the MUSIAL tool, which is used to build a storage object
		for genomic sequences. It writes the necessary parameters to a features file and a configuration file in
		the specified temporary directory. The configuration includes details about the reference genome, annotation,
		gene of interest, input VCF file, and thresholds for coverage and frequency.

	Arguments:
		reference (str): Path to the reference genome file.
		annotation (str): Path to the genome annotation file.
		gene (str): Name of the gene to be processed.
		vcf_input (str): Path to the input VCF file.
		tmpdir (str): Directory where temporary files will be written.
		min_coverage (int, optional): Minimum coverage threshold. Defaults to 3.
		min_frequency (float, optional): Minimum frequency threshold. Defaults to 0.8.

	Returns:
		None.

	Side Effects:
		Writes a features file (features.tsv) and a configuration file (config.json) to the specified temporary directory.
	"""
	with open( f"{tmpdir}/features.tsv", "w+" ) as features_file :
		features_file.write( f"gene\t{gene}\n" )
	config = {
		"minimalCoverage": min_coverage,
		"minimalFrequency": min_frequency,
		"storeFiltered": False,
		"skipSnpEff": True,
		"skipProteoformInference": False,
		"reference": reference,
		"annotation": annotation,
		"features": f"{tmpdir}/features.tsv",
		"output": f"{tmpdir}/storage.json",
		"vcfInput": [vcf_input]
	}
	with open( f"{tmpdir}/config.json", "w+" ) as config_file :
		json.dump( config, config_file )
		
def build_storage( musial_path : str, tmpdir : str ) :
	"""Builds the MUSIAL storage to extract genomic sequences from during a Nextstrain-TrepoGen gene workflow.

	Description:
		This function runs an external Java process to build a storage object for genomic sequences
		using the MUSIAL tool. It expects a configuration file (`config.json`) in the specified temporary directory.
		The process is executed with specific memory settings and the output is captured.

	Arguments:
		musial_path (str): Path to the MUSIAL JAR file.
		tmpdir (str): Path to the temporary directory containing 'config.json'.

	Returns:
		None.

	Side Effects:
		Writes the storage object with information about genomic sequences to a file in the specified temporary directory.
		
	Raises:
		ChildProcessError: If the subprocess returns a non-zero exit code, indicating an error during execution.
	"""
	process = subprocess.Popen(
			[
				"java",
				"-Xms1G",
				"-Xmx8G",
				"-jar",
				musial_path,
				"build",
				"-C",
				f"{tmpdir}/config.json"
			],
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
		)
	stdout, stderr = process.communicate()
	stdout = stdout.decode(encoding="utf-8")
	stderr = stderr.decode(encoding="utf-8")
	if process.returncode != 0:
		raise ChildProcessError(f"Error during execution of subprocess\nReturn code: {process.returncode}\n{stderr}")
	
def export_sequences( musial_path : str, gene : str, tmpdir : str ) :
	"""Exports gene sequences using MUSIAL during a Nextstrain-TrepoGen gene workflow.

	Description:
		This function runs an external Java process to extract nucleotide sequences for a specified gene,
		writing the output to a temporary directory. It expects the Java tool to produce exactly one `.fna`
		file, which is then renamed to `sequences.fasta` in the same directory.

	Arguments:
		musial_path (str): Path to the MUSIAL JAR file used for sequence extraction.
		gene (str): Name of the gene to extract sequences for.
		tmpdir (str): Path to a temporary directory for intermediate and output files.

	Returns:
		None.

	Side Effects:
		Writes a file named `sequences.fasta` in the specified temporary directory, containing the extracted gene sequences.

	Raises:
		ChildProcessError: If the Java subprocess reports an error (detected by "EXIT" in stderr).
		ValueError: If the number of `.fna` files produced is not exactly one.
	"""
	process = subprocess.Popen(
			[
				"java",
				"-Xms1G",
				"-Xmx8G",
				"-jar",
				musial_path,
				"sequence",
				"-c",
				"nt",
				"-k",
				"-r",
				"-F",
				gene,
				"-I",
				f"{tmpdir}/storage.json",
				"-o",
				f"{tmpdir}"
			],
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
		)
	stdout, stderr = process.communicate()
	stdout = stdout.decode(encoding="utf-8")
	stderr = stderr.decode(encoding="utf-8")
	if "EXIT" in stderr :
		raise ChildProcessError( f"Error during execution of subprocess\n{stderr}" )
	else :
		fna_files = glob(f"{tmpdir}/*.fna")
		if len(fna_files) != 1:
			raise ValueError(f"Expected exactly one .fna file, but found {len(fna_files)}: {fna_files}")
		os.rename(fna_files[0], f"{tmpdir}/sequences.fasta")

def get_insertions( gene_reference_sequence : str ) -> dict:
	"""Identifies and records insertions (gaps represented by '-') in a gene reference sequence.

	Description:
		This function scans a gene reference sequence for insertion gaps, which are represented by
		consecutive '-' characters. It returns a dictionary where the keys are the 1-based positions
		(relative to the ungapped sequence) where each insertion starts, and the values are the lengths
		of those insertions (number of consecutive '-' characters).

	Arguments:
		gene_reference_sequence (str): The gene reference sequence, where insertions are denoted by '-' characters.

	Returns:
		dict: A dictionary mapping the 1-based position (relative to the ungapped sequence) where each insertion starts
			  to the length of the insertion (number of consecutive '-' characters).

	Example:
		>>> get_insertions("ATG--C-T--A")
		{3: 2, 4: 1, 5: 2}
	"""
	insertions = { }
	position = 0
	start = 0
	length = 0
	cumulative_length = 0
	is_gap = False
	for char in gene_reference_sequence :
		if char == '-' :
			length += 1
			if not is_gap :
				start = position - cumulative_length
				is_gap = True
		elif is_gap :
			insertions[ start ] = length
			cumulative_length += length
			start = 0
			length = 0
			is_gap = False
		position += 1
	if is_gap :
		insertions[ start ] = length
	return insertions

def process_genome_annotation(annotation_file_path: str, gene_name: str, seq_record: Seq, region_names: set = set()) -> tuple[str, list, list, bool]:
	"""Processes a genome annotation file to extract and adjust features related to a specific gene.

	Description:
		This function reads a genome annotation file (typically in GFF format), locates the specified gene,
		and adjusts its features and associated biological regions to account for insertions and strand orientation.
		It also updates feature attributes and names for downstream analysis.

	Arguments:
		annotation_file_path (str): Path to the genome annotation file.
		gene_name (str): Name of the gene to extract and process.
		seq_record (Seq): Biopython Seq object representing the gene's sequence.
		region_names (set, optional): Set of region name prefixes to filter biological regions. Defaults to empty set.

	Returns:
		tuple:
			content (list of str): Adjusted annotation lines for the gene and its regions.
			names (list of str): List of processed biological region names.
			biological_regions (list of tuple): List of tuples containing (region name, start, end) for each biological region.
			is_reverse (bool): True if the gene is on the reverse strand, False otherwise.
	Notes:
		- The function assumes the presence of a global variable `biological_region_name_map` for name replacements.
		- The function relies on an external `get_insertions` function to identify insertion-induced gaps.
		- The function modifies the strand orientation and coordinates of features as needed.
	"""
	
	with open( annotation_file_path, "r" ) as genome_annotation_file :
		lines = genome_annotation_file.readlines()
	content = [] # Initialize content for the adjusted annotation file.
	start = None # Start position of the gene feature.
	end = None # End position of the gene feature.
	locus_tag = None # Locus tag of the gene feature.
	is_reverse = False # Flag to indicate if the gene is on the reverse strand.
	insertions = {} # Dictionary to store insertion-induced gaps in the reference sequence.
	cumulative_insertion_length = None # Function to calculate cumulative insertion length at a given position.
	names = [] # List to store names of gene sub-features.
	biological_regions = [] # List to store biological regions associated with the gene feature.
	global biological_region_name_map
	
	for line in lines : # Iterate through lines in the annotation file until specific gene feature is found.
		fields = []
		if line.startswith( '#' ) : # Line is a comment or header line.
			if line.startswith( '##sequence-region' ) : # Adjust the sequence region header line to the length of the gene's sequence.
				fields = line.strip().split(' ')
				content.append( ' '.join( [ fields[0], fields[1], str(1), str(len(seq_record)) ] ) )
			else :
				content.append( line.strip() )
			continue
		else :
			fields = line.strip().split( '\t' )

		if fields[2] == "gene" and f"gene={gene_name}" in line : # Line contains the full gene feature.
			start = int(fields[3])
			end = int(fields[4])
			attributes = dict( [ _.split( "=" ) for _ in fields[8].split(";") ] )
			locus_tag = attributes.get( 'locus_tag', None )
			is_reverse = fields[6] == '-'
			if is_reverse :
				seq_record = seq_record.reverse_complement( )
			insertions = get_insertions( str( seq_record ) ) # Get insertion-induced gaps in the reference sequence.
			cumulative_insertion_length = lambda p: sum( [ l for s, l in insertions.items() if s < p ] )

	for line in filter(lambda l : not l.startswith("#") and l.strip().split('\t')[2] == 'biological_region' and f"locus_tag={locus_tag}" in l, lines) :
		fields = line.strip().split( '\t' )
		# Adjust start and end for insertion induced gaps and strand orientation.
		if is_reverse :
			s = end - int(fields[4]) + 1
			e = end - int(fields[3]) + 1
			fields[6] = '+' # Reverse strand orientation of entry.
		else :
			s = int(fields[3]) - start + 1
			e = int(fields[4]) - start + 1
		s += cumulative_insertion_length( s )
		e += cumulative_insertion_length( e )
		fields[3] = str(s)
		fields[4] = str(e)

		# Adjust feature type to gene.
		fields[2] = 'gene'
		
		# Adjust feature attributes.
		attributes = dict( [ _.split( "=" ) for _ in fields[8].split(";") ] )
		attributes["locus_tag"] = attributes[ "ID" ]
		attributes["ID"] = 'gene-' + attributes["ID"]
		del attributes["Parent"]
		attributes["gbkey"] = 'gene'
		name = attributes['gene'].split('-')[0] # Only keep the biological region prefix as the name.
		if any( [ name.startswith( r ) for r in region_names ] ) :
			biological_regions.append( ( name, fields[3], fields[4] ) )
		
		# Replace biological region names with more descriptive names.
		for k, v in biological_region_name_map.items( ) :
			name = name.replace( k, v )
		names.append( name )
		
		# Replace the gene name in the attributes and re-convert attributes to a string.
		attributes[ 'gene' ] = name
		attributes[ 'Name' ] = name
		fields[8] = ';'.join( [ str(k) + '=' + str(v) for k, v in attributes.items( ) ] )
		
		# Write the adjusted line to the content.
		content.append( '\t'.join( fields ) )

	return content, names, biological_regions, is_reverse

def translate_biological_regions(sample_records: list, region_type: str, biological_regions: list) -> tuple[dict, dict, dict, dict]:
	"""Translates specified biological regions from a list of sample records and organizes the results.

	Description:
		This function processes a list of SeqRecord objects, extracting and translating sequences
		for specified biological regions that match the given region_type. It returns several mappings:
		- A mapping of sample names to their translated sequences for each region.
		- A mapping of sample names to their concatenated region type indices.
		- A count of how many samples have each type.
		- A mapping of region names to their unique translated sequence indices.

	Arguments:
		sample_records (list): A list of SeqRecord objects, each representing a biological sample with an 'id' and 'seq'.
		region_type (str): The prefix of the region names to process (e.g., "ecl", "sp").
		biological_regions (list): A list of regions, where each region is a tuple or list containing:
			- region[0]: The region name (str).
			- region[1]: The start position (int or str convertible to int, 1-based inclusive).
			- region[2]: The end position (int or str convertible to int, 1-based inclusive).

	Returns:
		tuple:
			samples_per_region_sequence (dict): 
				Mapping of sample name to a dict of region name to translated amino acid sequence.
			samples_type (dict): 
				Mapping of sample name to a string representing the concatenated region type indices.
			types_count (dict): 
				Mapping of region type string to the number of samples with that type.
			regions_per_sequence_index (dict): 
				Mapping of region name to a dict of translated sequence to its unique index.
	
	Notes:
		- Only regions whose names start with the specified region_type are processed.
		- Sequences are padded with 'N' to ensure their length is a multiple of 3 before translation.
		- Translation is performed using Bio.Seq and Bio.SeqRecord.
		- The function assumes that sample_records are compatible with Bio.SeqRecord and that SeqIO and Seq are imported.
	"""
	samples_per_region_sequence = {}
	samples_type = {}
	types_count = {}
	regions_per_sequence_index = {}

	for sample_record in sample_records :
		sample_name = str( sample_record.id )
		sample_sequence = str( sample_record.seq )
		samples_per_region_sequence.setdefault( sample_name, {} )
		
		for region in biological_regions :
			if region[0].startswith(region_type) : # Only process regions of the specified type.
				region_name = region[0]
				_sequence = sample_sequence[ int(region[1]) - 1 : int(region[2]) ].replace( '-', '' ) 
				if len(_sequence) % 3 != 0: # Ensure the sequence length is a multiple of 3 for translation.
					_sequence += 'N' * (3 - len(_sequence) % 3)
				_translated_sequence = str( SeqIO.SeqRecord( Seq( _sequence ) ).translate().seq )

				regions_per_sequence_index.setdefault( region_name, {} )
				if not _translated_sequence in regions_per_sequence_index[ region_name ] : # If the translated sequence is not already indexed, add it.
					regions_per_sequence_index[ region_name ][ _translated_sequence ] = len( regions_per_sequence_index[ region_name ] )
				samples_per_region_sequence[ sample_name ][ region_name ] = _translated_sequence

		_type = '.'.join( [ str( regions_per_sequence_index[r][s] ) for r, s in samples_per_region_sequence[ sample_name ].items( ) ] )
		samples_type[ sample_name ] = _type
		types_count[ _type ] = types_count.get( _type, 0 ) + 1

	return samples_per_region_sequence, samples_type, types_count, regions_per_sequence_index

def main():
	"""Extract gene specific data for a Nextstrain-TrepoGen gene dataset.

	Description:
		Orchestrates the workflow for extracting gene specific sequences and annotations,
		processing gene sub-features, and generating output files required for further
		analysis and visualization. The workflow:
		- Builds a MUSIAL alignment for the target gene (including its config file).
		- Parses the genome annotation to extract and type gene sub-features (e.g., topologies,
		active sites).
		- Exports the reformatted reference gene sequence, a gap-adjusted GFF3 with updated
		coordinates (handling reverse complements as needed), and a list of sub-feature names.
		- Assigns sub-feature types to each sample, producing per-sample metadata.
		- Generates an Auspice config for Nextstrain (custom color schemes).

	Arguments:
		None. Arguments are parsed internally via `parse_args()`.
	
	Returns:
		None.

	Side Effects:
		Writes several output files to the specified paths, including:
		- Reference gene sequence (FASTA).
		- Updated annotation (GFF3).
		- Gene names (TXT).
		- Per sample region type metadata (CSV).
		- Sample sequences (FASTA).
		- Auspice configuration (JSON).
	
	Raises:
		ValueError: If the reference sequence identifier cannot be parsed from the reference FASTA file.

	Note:
		This function assumes the existence of several helper functions and external dependencies,
		such as MUSIAL, Biopython's SeqIO, and pandas.
	"""
	
	args = parse_args()

	with tempfile.TemporaryDirectory( ) as tmpdir :
		# Write configuration files for MUSIAL.
		write_config( args.reference, args.annotation, args.gene, args.vcf, tmpdir )

		# Build storage for genomic sequences using the specified MUSIAL path and working directory.
		build_storage( args.musial, tmpdir )
	
		# Export sequences for the specified gene feature using MUSIAL.
		export_sequences( args.musial, args.gene, tmpdir )

		# Parse the reference (genome) sequence identifier to be reused in the reference (gene) sequence file.
		reference_identifier = None
		for record in SeqIO.parse( args.reference, "fasta" ) :
			reference_identifier = str( record.id )
			break
		if reference_identifier is None :
			raise ValueError( f"Could not parse reference sequence identifier from {args.reference}" )

		# Parse per sample sequence records (the first record is expected to be the reference sequence of the specified gene).
		samples_sequence_records = SeqIO.to_dict( SeqIO.parse( f"{tmpdir}/sequences.fasta", "fasta" ) )
		reference_seq_record = samples_sequence_records[ "reference" ].seq

		# Process genome annotation; extract sub-features of the specified gene and biological regions for typing.
		gene_annotation_content, sub_feature_names, biological_regions, is_reverse = process_genome_annotation(args.annotation, args.gene, reference_seq_record, set( [ 'ecl' ] ) )
		if is_reverse :
			for sample_sequence_record in samples_sequence_records.values( ):
				sample_sequence_record.seq = sample_sequence_record.seq.reverse_complement( )

		# Write the extracted reference sequence to the output file.
		with open( f"{args.output_reference}", "w+" ) as gene_reference_sequence_file :
			gene_reference_sequence_file.write( f">{reference_identifier} [gene={args.gene}]\n" )
			gene_reference_sequence_file.write( "\n".join( wrap( str( samples_sequence_records["reference"].seq ), 60, break_on_hyphens=False ) ) + "\n" )
		# Write the updated annotation content and gene's sub-feature names.
		with open( f"{args.output_annotation}", "w+" ) as gff_output_file :
			gff_output_file.write( "\n".join( gene_annotation_content ) + "\n" )
		with open( f"{args.output_genes}", "w+" ) as genes_output_file :
			genes_output_file.write( "\n".join( sub_feature_names ) + "\n" )

		# Analyze the sample sequences and biological regions to create metadata for specified region's types.
		samples_type_records = {}
		auspice_config = {"colorings": []}
		global biological_region_name_map
		for region in [ 'ecl' ] : # Currently, only ECL regions are processed.
			samples_per_region_sequence, samples_type, types_count, regions_per_sequence_index = translate_biological_regions( samples_sequence_records.values(), region, biological_regions )
			
			n = len( samples_type ) # Number of samples.
			_types = set( [ ] ) # Set of actually exported types; filtered types are not included.
			for _name, _type in samples_type.items() :
				samples_type_records.setdefault( _name, {"sra": _name} )
				if _type != samples_type["reference"] and types_count[ _type ] / n < 0.01 :
					_type = 'other' # Replace rare ECL types with 'other'.
				else :
					_types.add( _type )
				samples_type_records[ _name ][ f"{region}_type" ] = _type
				for _region, _sequence in samples_per_region_sequence[ _name ].items( ) :
					samples_type_records[ _name ][ _region ] = _sequence

			# Add custom region type coloring.
			_types = sorted( list( _types ), key = lambda t : types_count.get( t, 0 ), reverse = True )
			_colors = [ colors.rgb2hex(c) for c in colormaps['managua'](linspace(0, 1, len(_types))) ]
			scale = [ ]
			for idx, _type in enumerate( _types ) :
				if _type == 'other' :
					scale.append( [ _type, "#aaaaaa" ] )
				else :
					scale.append( [ _type, _colors[idx] ] )
			auspice_config[ "colorings" ].append(
				{
					"key": f"{region}_type",
					"title": f"{biological_region_name_map[region].strip().replace('-',' ')} type",
					"type": "categorical",
					"scale": scale,
				}
			)

			# Add custom region sequence coloring.
			for _region, _ in regions_per_sequence_index.items() :
				scale = [ ]
				legend = [ ]
				_colors = [ colors.rgb2hex(c) for c in colormaps['managua'](linspace(0, 1, len(_))) ]
				for _sequence, _index in _.items() :
					scale.append( [ _sequence, _colors[_index] ] )
					legend.append( { "value": _sequence, "display": _index } )
				auspice_config[ "colorings" ].append(
					{
						"key" : _region,
						"title" : f"{_region.replace(region,biological_region_name_map[region]).strip().replace('-',' ')}",
						"type" : "categorical",
						"scale": scale,
						"legend": legend,
					}
				)
					
		# If the defined root identifier is present in the records, remove it; otherwise, set it as the identifier for the reference sequence.
		if args.root is not None :
			if args.root in samples_sequence_records :
				del samples_type_records[ "reference" ]
				del samples_sequence_records[ "reference" ]
			else :
				samples_type_records[ "reference" ][ "sra" ] = args.root
				samples_sequence_records[ "reference" ].id = args.root

		with open( f"{args.output_types}", "w+" ) as types_output_file :
			df = DataFrame.from_records( list( samples_type_records.values() ) )
			df.to_csv( types_output_file, index=False )

		# Write the adjusted auspice configuration file.
		with open( args.output_configuration, "w+" ) as config_file :
			json.dump( auspice_config, config_file, indent=2 )

		# Write the updated records back to the sequences FASTA file.
		SeqIO.write( samples_sequence_records.values(), f"{args.output_sequences}", "fasta" )

if __name__ == "__main__":
	main()
