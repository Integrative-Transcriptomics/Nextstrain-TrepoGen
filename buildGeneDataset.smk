configfile: "config/config.yaml"

# sources: List of data sources to be used in the workflow. Each entry represents a dataset of source/data.
sources = [
	"tpass308"
]
# subsets: List of variant subsets to process. Possible values include "snv" (single nucleotide variants) and "snv-indel" (single nucleotide variants and insertions/deletions).
subsets = [
	"snv",
	"snv-indel"
]
# genes: List of gene identifiers to include in the analysis. Each entry corresponds to a gene of interest.
genes = [
	"Tp0865"
]

# data_source_files: List of files that are part of the data source. These files are expected to be present in the source/data/{source} directory.
data_source_files = [
	"sequence.fasta", # Reference genome sequence in FASTA format.
	"annotation.gff3", # Annotation file in GFF3 format.
	"meta.csv", # Metadata file in CSV format containing information about the single samples.
	"meta_colors.tsv", # Metadata file in TSV format containing color information.
	"configuration.json", # Auspice configuration file in JSON format providing settings for the dataset.
	"description.md" # Auspice description file in Markdown format (MD) providing details about the dataset.
]

# build_files: List of files that are part of the build process. These files are expected to be present in the builds/{source}_{subset}_{gene} directory.
build_files = [
	"exclude.txt", # Exclusion file listing samples to be excluded from the analysis.
	"mask.bed", # BED file specifying regions to be masked in the genome.
	"configuration.json", # Auspice configuration file in JSON format for the build.
	"description.md" # Auspice description file in Markdown format (MD) for the build.
]

# work_files: List of files that are generated during the workflow. These files will be generated in the builds/{source}_{subset}_{gene}/.work directory.
work_files = [
	"reference.fasta", # Reference gene sequence in FASTA format.
	"annotation.gff3", # Annotation file in GFF3 format for the gene.
	"names.txt", # Text file containing names of features in the gene.
	"sequences.fasta", # Sample sequences in FASTA format for the gene.
	"types.csv", # Types file in CSV format containing information about selected feature types, e.g. extracellular domains.
	"types_configuration.json", # Auspice configuration file in JSON format for the selected feature types.
	"meta.csv", # Metadata file in CSV format containing information about the samples.
	"index.tsv", # Index file in TSV format for the sequences.
	"filtered.fasta", # Filtered sequences in FASTA format after applying the filter rule.
	"initial.nwk", # Initial tree in Newick format (NWK) before refinement.
	"tree.nwk", # Refined tree in Newick format (NWK) after refinement.
	"auspice.json", # Auspice JSON file containing the final dataset for visualization.
	"branch_lengths.json", # Branch lengths in JSON format for the tree.
	"nucleotide_mutations.json", # Nucleotide mutations in JSON format for the tree.
	"traits.json", # Traits in JSON format for the tree.
	"amino_acid_mutations.json", # Amino acid mutations in JSON format for the tree.
	"configuration.json", # Final Auspice configuration file in JSON format for the dataset.
	"description.md", # Final Auspice description file in Markdown format (MD) for the dataset.
]

# Defines the directory structure and I/O files for each source, subset, and gene.
rule all:
	input:
		expand("source/data/{source}/{data_source_file}", source=sources, data_source_file=data_source_files),
		expand("source/data/{source}/variants.{subset}.vcf", source=sources, subset=subsets),
		expand("builds/{source}_{subset}_{gene}/{build_file}", source=sources, gene=genes, subset=subsets, build_file=build_files),
		expand("builds/{source}_{subset}_{gene}/.work/{work_file}", source=sources, gene=genes, subset=subsets, work_file=work_files),
		expand("datasets/{source}_{subset}_{gene}.json", source=sources, gene=genes, subset=subsets)

# Prepares the gene build by processing the reference genome, annotation, and variants into gene specific data using MUSIAL. 
rule prepare:
	input:
		reference="source/data/{source}/sequence.fasta",
		annotation="source/data/{source}/annotation.gff3",
		variants="source/data/{source}/variants.{subset}.vcf",
	output:
		reference="builds/{source}_{subset}_{gene}/.work/reference.fasta",
		annotation="builds/{source}_{subset}_{gene}/.work/annotation.gff3",
		features="builds/{source}_{subset}_{gene}/.work/names.txt",
		sequences="builds/{source}_{subset}_{gene}/.work/sequences.fasta",
		types="builds/{source}_{subset}_{gene}/.work/types.csv",
		configuration="builds/{source}_{subset}_{gene}/.work/types_configuration.json",
	params:
		root=lambda wc: config.get(wc.source, {}).get("reference_sample", ""),
	shell:
		"""
		python scripts/prepare_gene_build.py \
			-ir {input.reference} \
			-ia {input.annotation} \
			-iv {input.variants} \
			-ig {wildcards.gene} \
			-or {output.reference} \
			-oa {output.annotation} \
			-og {output.features} \
			-os {output.sequences} \
			-ot {output.types} \
			-oc {output.configuration} \
			-m /mnt/c/Users/siha/Desktop/Tools/MUSIAL-v2.4.1.jar \
			-q {params.root}
		"""

# Merges the description files from the source and build into a single description file for the dataset.
rule merge_description:
	input:
		build_description="builds/{source}_{subset}_{gene}/description.md",
		source_description="source/data/{source}/description.md",
	output:
		"builds/{source}_{subset}_{gene}/.work/description.md",
	shell:
		"""
		python scripts/concat_files.py -o {output} {input}
		"""

# Merges the configuration files from the source and build into a single configuration file for the dataset.
rule merge_configuration:
	input:
		source_configuration="source/data/{source}/configuration.json",
		build_configuration="builds/{source}_{subset}_{gene}/configuration.json",
		types_configuration="builds/{source}_{subset}_{gene}/.work/types_configuration.json",
	output:
		"builds/{source}_{subset}_{gene}/.work/configuration.json",
	shell:
		"""
		python scripts/concat_files.py -o {output} {input}
		"""

# Merges the metadata files from the source and build into a single metadata file for the dataset.
rule merge_metadata:
	input:
		metadata_source="source/data/{source}/meta.csv",
		metadata_types="builds/{source}_{subset}_{gene}/.work/types.csv",
	output:
		"builds/{source}_{subset}_{gene}/.work/meta.csv"
	params:
		metadata_id=lambda wc: config.get(wc.source, {}).get("meta_identifier", "name strain id"),
	shell:
		"""
		augur merge \
			--metadata source={input.metadata_source} types={input.metadata_types} \
			--metadata-id-columns {params.metadata_id} \
			--output-metadata {output}
		"""

# Merges the color files from the source and geo into a single color file for the dataset.
rule merge_color:
	input:
		source_colors="source/data/{source}/meta_colors.tsv",
		geo_colors="source/geo/color.tsv",
	output:
		"builds/{source}_{subset}_{gene}/.work/colors.tsv",
	shell:
		"""
		python scripts/concat_files.py -o {output} {input}
		"""

rule index:
	input:
		rules.prepare.output.sequences,
	output:
		"builds/{source}_{subset}_{gene}/.work/index.tsv",
	shell:
		"""
		augur index \
		  --sequences {input} \
		  --output {output}
		"""

# Filters the sequences based on the metadata and excludes specified samples. If no filter is defined, it simply passes the sequences through.
# TODO: This may be applied to the variants before the prepare step.
rule filter:
	input:
		sequences=rules.prepare.output.sequences,
		sequence_index=rules.index.output,
		metadata="builds/{source}_{subset}_{gene}/.work/meta.csv",
		exclude="builds/{source}_{subset}_{gene}/exclude.txt",
	output:
		"builds/{source}_{subset}_{gene}/.work/filtered.fasta",
	params:
		query=lambda wc: config.get(wc.source, {}).get("filter.query", None),
		metadata_id=lambda wc: config.get(wc.source, {}).get("meta_identifier", "name strain id"),
	run:
		if params.query:
			shell(
				"""
				augur filter --sequences {input.sequences} \
					--sequence-index {input.sequence_index} \
					--metadata {input.metadata} \
					--metadata-id-columns {params.metadata_id} \
					--exclude {input.exclude} \
					--output-sequences {output} \
					{params.query} \
				"""
			)
		else:
			# If no filter is defined, just pass the file through.
			shell("cp {input.sequences} {output}")

# Builds the phylogenetic tree from the filtered alignment.
rule tree:
	input:
		alignment=rules.filter.output,
		exclude="builds/{source}_{subset}_{gene}/mask.bed",
	output:
		"builds/{source}_{subset}_{gene}/.work/initial.nwk",
	params:
		method=config.get( "gene", {} ).get( "tree.method", "--method iqtree" ),
	shell:
		"""
		augur tree --alignment {input.alignment} \
			--exclude-sites {input.exclude} \
			--output {output} \
			{params.method}
		"""

# Refines the phylogenetic tree; performs date inference and branch length estimation.
rule refine:
	input:
		alignment=rules.filter.output,
		metadata="builds/{source}_{subset}_{gene}/.work/meta.csv",
		tree="builds/{source}_{subset}_{gene}/.work/initial.nwk",
	output:
		tree="builds/{source}_{subset}_{gene}/.work/tree.nwk",
		branch_lengths="builds/{source}_{subset}_{gene}/.work/branch_lengths.json",
	params:
		metadata_id=lambda wc: config.get(wc.source, {}).get("meta_identifier", "name strain id"),
		clock_rate=lambda wc: config.get(wc.gene, {}).get("refine.clock_rate", ""),
		iterations=config.get("refine.iterations", 1),
		precision=config.get("refine.precision", 1),
		root=lambda wc: config.get(wc.source, {}).get("refine.root", ""),
		year_bounds=lambda wc: config.get(wc.source, {}).get("refine.year_bounds", ""),
		seed=config.get("seed", 1),
	shell:
		"""
		augur refine --tree {input.tree} \
			--alignment {input.alignment} \
			--metadata {input.metadata} \
			--metadata-id-columns {params.metadata_id} \
			--timetree \
			--max-iter {params.iterations} \
			--covariance \
			--keep-polytomies \
			--precision {params.precision} \
			--coalescent opt \
			--date-format %Y-%m-%d \
			--date-inference marginal \
			--date-confidence \
			--branch-length-inference marginal \
			--seed {params.seed} \
			--output-tree {output.tree} \
			--output-node-data {output.branch_lengths} \
			{params.root} \
			{params.clock_rate} \
			{params.year_bounds}
		"""

# Computes the ancestral sequences for the phylogenetic tree, inferring nucleotide mutations.
rule ancestral:
	input:
		alignment=rules.filter.output,
		reference=rules.prepare.output.reference,
		tree="builds/{source}_{subset}_{gene}/.work/tree.nwk",
	output:
		"builds/{source}_{subset}_{gene}/.work/nucleotide_mutations.json",
	params:
		seed=config.get("seed", 1),
	shell:
		"""
		augur ancestral --tree {input.tree} \
			--alignment {input.alignment} \
			--root-sequence {input.reference} \
			--inference marginal \
			--keep-ambiguous \
			--keep-overhangs \
			--seed {params.seed} \
			--output-node-data {output}
		"""

# Computes the ancestral traits for the phylogenetic tree, inferring traits from the metadata.
rule traits:
	input:
		tree="builds/{source}_{subset}_{gene}/.work/tree.nwk",
		metadata="builds/{source}_{subset}_{gene}/.work/meta.csv",
	output:
		"builds/{source}_{subset}_{gene}/.work/traits.json",
	params:
		metadata_id=lambda wc: config.get(wc.source, {}).get("meta_identifier", "name strain id"),
		columns=lambda wc: config.get(wc.source, {}).get("traits.columns", ""),
	shell:
		"""
		augur traits --tree {input.tree} \
			--metadata {input.metadata} \
			--metadata-id-columns {params.metadata_id} \
			--columns {params.columns} \
			--output-node-data {output}
		"""

# Translates nucleotide mutations of specified gene features into amino acid mutations.
rule translate:
	input:
		tree="builds/{source}_{subset}_{gene}/.work/tree.nwk",
		sequences=rules.ancestral.output,
		annotation=rules.prepare.output.annotation,
	output:
		"builds/{source}_{subset}_{gene}/.work/amino_acid_mutations.json",
	params:
		genes=rules.prepare.output.features,
	shell:
		"""
		augur translate \
			--tree {input.tree} \
			--ancestral-sequences {input.sequences} \
			--reference-sequence {input.annotation} \
			--output-node-data {output} \
			--genes {params.genes}
		"""

# Exports the refined tree and associated data into a format suitable for visualization in Auspice.
rule export:
	input:
		tree=rules.refine.output.tree,
		metadata="builds/{source}_{subset}_{gene}/.work/meta.csv",
		colors=rules.merge_color.output,
		configuration="builds/{source}_{subset}_{gene}/.work/configuration.json",
		description="builds/{source}_{subset}_{gene}/.work/description.md",
		coordinates="source/geo/loc.tsv",
		branch_lengths="builds/{source}_{subset}_{gene}/.work/branch_lengths.json",
		traits="builds/{source}_{subset}_{gene}/.work/traits.json",
		nucleotide_mutations="builds/{source}_{subset}_{gene}/.work/nucleotide_mutations.json",
		amino_acid_mutations="builds/{source}_{subset}_{gene}/.work/amino_acid_mutations.json",
	output:
		"builds/{source}_{subset}_{gene}/.work/auspice.json",
	params:
		metadata_id=lambda wc: config.get(wc.source, {}).get("meta_identifier", "name strain id"),
	shell:
		"""
		augur export v2 \
			--tree {input.tree} \
			--metadata {input.metadata} \
			--metadata-id-columns {params.metadata_id} \
			--node-data {input.branch_lengths} {input.traits} {input.nucleotide_mutations} {input.amino_acid_mutations} \
			--auspice-config {input.configuration} \
			--description {input.description} \
			--colors {input.colors} \
			--lat-longs {input.coordinates} \
			--output {output} \
			--include-root-sequence
		"""

# Reprocesses the Auspice JSON file; currently this only applies a color adjustment to the Auspice dataset.
rule reprocess_auspice:
	input:
		"builds/{source}_{subset}_{gene}/.work/auspice.json",
	output:
		"datasets/{source}_{subset}_{gene}.json",
	shell:
		"""
		python scripts/reprocess_auspice.py \
			-i {input} \
			-o {output}
		"""
