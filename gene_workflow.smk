configfile: "config/global.yaml"
configfile: "config/data.yaml"
configfile: "config/genes.yaml"

# genes: List of genes to use in the workflow; see config/config.yaml.
genes = list(config.get("genes").keys())

# sources: List of source datasets to use in the workflow. Expected to match one in source/data/<VAL>.
sources = ["TPASS-308"]

# subsets: List of variant subsets to use in the workflow. Expected to match one in source/data/*/variants/<VAL>.vcf
subsets = ["snv"]

# source_files: List of files that are part of the source data. These files are expected to be present in the source/data/{source} directory.
source_files = [
	"sequence.fasta", # Reference genome sequence in FASTA format.
	"annotation.gff3", # Annotation file in GFF3 format.
	"meta.csv", # Metadata file in CSV format containing information about the single samples.
	"meta_colors.tsv", # Metadata file in TSV format containing color information.
	"auspice_configuration.json", # Auspice configuration file in JSON format providing settings for the dataset.
	"auspice_description.md" # Auspice description file in Markdown format (MD) providing details about the dataset.
]

# work_files: List of files that are generated during the workflow. These files will be generated in the .work/{source}_{subset}_{gene}/ directory.
work_files = [
	"reference.fasta", # Reference gene sequence in FASTA format.
	"annotation.gff3", # Annotation file in GFF3 format for the gene.
	"gene_features.txt", # Text file containing names of features in the gene.
	"sequences.fasta", # Sample sequences in FASTA format for the gene.
	"gene_meta.csv", # Types file in CSV format containing information about selected feature types, e.g. extracellular domains.
	"gene_auspice_configuration.json", # Auspice configuration file in JSON format for the selected feature types.
	"meta.csv", # Metadata file in CSV format containing information about the samples.
	"initial.nwk", # Initial tree in Newick format (NWK) before refinement.
	"tree.nwk", # Refined tree in Newick format (NWK) after refinement.
	"auspice.json", # Auspice JSON file containing the final dataset for visualization.
	"branch_lengths.json", # Branch lengths in JSON format for the tree.
	"nucleotide_mutations.json", # Nucleotide mutations in JSON format for the tree.
	"traits.json", # Traits in JSON format for the tree.
	"amino_acid_mutations.json", # Amino acid mutations in JSON format for the tree.
]

# Defines the directory structure and I/O files for each source, subset, and gene.
rule all:
	input:
		expand("source/data/{source}/{source_file}", source=sources, source_file=source_files),
		expand("source/data/{source}/variants/{subset}.vcf", source=sources, subset=subsets),
		expand("source/data/{source}/variants/{subset}.tsv", source=sources, subset=subsets),
		expand("source/data/{source}/variants/{subset}.filtered.vcf", source=sources, subset=subsets),
		expand(".work/{source}_{subset}_{gene}/{work_file}", source=sources, gene=genes, subset=subsets, work_file=work_files),
		"source/geo/color.tsv", # Independent file for geo colors.
		"source/geo/loc.tsv" # Independent file for geo coordinates.

# Generates an index of the input variants.
rule index:
	input:
		"source/data/{source}/variants/{subset}.vcf",
	output:
		"source/data/{source}/variants/{subset}.tsv",
	shell:
		"""
		augur index \
		  --sequences {input} \
		  --output {output}
		"""

# Filters variants (optional) based on metadata and excludes specified samples. If no filter is defined, it simply passes the data through.
rule filter:
	input:
		variants="source/data/{source}/variants/{subset}.vcf",
		index="source/data/{source}/variants/{subset}.tsv",
		metadata="source/data/{source}/meta.csv",
	output:
		"source/data/{source}/variants/{subset}.filtered.vcf",
	params:
		metadata_id=lambda wc: config.get(wc.source).get("meta_identifier", "name strain id"),
		query_cl=lambda wc: config.get(wc.source).get("filter.query_cl", ""),  # Optional query command line argument in config to filter variants.
		exclude_cl=lambda wc: config.get(wc.source).get("filter.exclude_cl", ""),  # Optional exclude comman line argument to exclude samples.
	run:
		if bool(params.query_cl) or bool(params.exclude_cl):
			shell(
				"""
				augur filter --sequences {input.variants} \
					--sequence-index {input.index} \
					--metadata {input.metadata} \
					--metadata-id-columns {params.metadata_id} \
					--output-sequences {output} \
					{params.query_cl} \
					{params.exclude_cl}
				"""
			)
		else:
			# If no filter is defined, just pass the file through.
			shell("cp {input.variants} {output}")

# Prepares a gene build by processing the reference genome, annotation, and variants into gene specific data.
rule prepare:
	input:
		reference="source/data/{source}/sequence.fasta",
		annotation="source/data/{source}/annotation.gff3",
		variants=rules.filter.output,
	output:
		reference=".work/{source}_{subset}_{gene}/reference.fasta",
		annotation=".work/{source}_{subset}_{gene}/annotation.gff3",
		features=".work/{source}_{subset}_{gene}/gene_features.txt",
		sequences=".work/{source}_{subset}_{gene}/sequences.fasta",
		gene_meta=".work/{source}_{subset}_{gene}/gene_meta.csv",
		gene_auspice_configuration=".work/{source}_{subset}_{gene}/gene_auspice_configuration.json",
	params:
		musial=config.get("musial"),
		root=lambda wc: config.get(wc.source).get("reference_sample", ""),
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
			-ot {output.gene_meta} \
			-oc {output.gene_auspice_configuration} \
			-m {params.musial} \
			-q {params.root}
		"""

# Merges metadata files of the source data and the gene build.
rule metadata:
	input:
		metadata_source="source/data/{source}/meta.csv",
		metadata_gene=rules.prepare.output.gene_meta,
	output:
		".work/{source}_{subset}_{gene}/meta.csv"
	params:
		metadata_id=lambda wc: config.get(wc.source).get("meta_identifier", "name strain id"),
	shell:
		"""
		augur merge \
			--metadata msource={input.metadata_source} mgene={input.metadata_gene} \
			--metadata-id-columns {params.metadata_id} \
			--output-metadata {output}
		"""

# Concatenates metadata color definition files of the source data and the geo information.
rule colors:
	input:
		source_colors="source/data/{source}/meta_colors.tsv",
		geo_colors="source/geo/color.tsv",
	output:
		".work/{source}_{subset}_{gene}/colors.tsv",
	shell:
		"""
		cat {input.source_colors} {input.geo_colors} >> {output}
		"""

# Builds the phylogenetic tree from the sequence alignment.
rule tree:
	input:
		alignment=rules.prepare.output.sequences,
	output:
		".work/{source}_{subset}_{gene}/initial.nwk",
	params:
		method_cl=lambda wc: config.get("genes").get(wc.gene).get("tree.method_cl", "--method iqtree"),
		exclude_cl=lambda wc: config.get("genes").get(wc.gene).get("tree.exclude_cl", ""),
	shell:
		"""
		augur tree --alignment {input.alignment} \
			--output {output} \
			{params.method_cl} \
			{params.exclude_cl}
		"""

# Refines the phylogenetic tree; performs date inference and branch length estimation.
rule refine:
	input:
		alignment=rules.prepare.output.sequences,
		metadata=rules.metadata.output,
		tree=rules.tree.output,
	output:
		tree=".work/{source}_{subset}_{gene}/tree.nwk",
		branch_lengths=".work/{source}_{subset}_{gene}/branch_lengths.json",
	params:
		seed=config.get("seed", 1),
		iterations=config.get("refine.iterations", 1),
		precision=config.get("refine.precision", 1),
		metadata_id=lambda wc: config.get(wc.source).get("meta_identifier", "name strain id"),
		clock_rate_cl=lambda wc: config.get("genes").get(wc.gene).get("refine.clock_rate_cl", ""),
		root_cl=lambda wc: config.get(wc.source).get("refine.root_cl", ""),
		year_bounds_cl=lambda wc: config.get(wc.source).get("refine.year_bounds_cl", ""),
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
			{params.root_cl} \
			{params.clock_rate_cl} \
			{params.year_bounds_cl}
		"""

# Computes the ancestral sequences for the phylogenetic tree, inferring nucleotide mutations.
rule ancestral:
	input:
		alignment=rules.prepare.output.sequences,
		reference=rules.prepare.output.reference,
		tree=rules.refine.output.tree,
	output:
		".work/{source}_{subset}_{gene}/nucleotide_mutations.json",
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
		tree=rules.refine.output.tree,
		metadata=rules.metadata.output,
	output:
		".work/{source}_{subset}_{gene}/traits.json",
	params:
		metadata_id=lambda wc: config.get(wc.source).get("meta_identifier", "name strain id"),
		columns=lambda wc: config.get(wc.source).get("traits.columns", ""),
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
		tree=rules.refine.output.tree,
		sequences=rules.ancestral.output,
		annotation=rules.prepare.output.annotation,
	output:
		".work/{source}_{subset}_{gene}/amino_acid_mutations.json",
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
		metadata=rules.metadata_metadata.output,
		colors=rules.colors.output,
		source_auspice_configuration="source/data/{source}/auspice_configuration.json",
		gene_auspice_configuration=rules.prepare.output.gene_auspice_configuration,

		description="source/data/{source}/auspice_description.md",

		coordinates="source/geo/loc.tsv",
		branch_lengths=rules.refine.output.branch_lengths,
		traits=rules.traits.output,
		nucleotide_mutations=rules.ancestral.output,
		amino_acid_mutations=rules.translate.output,
	output:
		".work/{source}_{subset}_{gene}/auspice.json",
	params:
		metadata_id=lambda wc: config.get(wc.source).get("meta_identifier", "name strain id"),
		title=lambda wc: f"'{config.get('export.title', 'TrepoGen')} ({wc.gene})'",
		maintainers=config.get("export.maintainers", ""),
		build_url=config.get("export.build_url", ""),
	shell:
		"""
		augur export v2 \
			--tree {input.tree} \
			--metadata {input.metadata} \
			--metadata-id-columns {params.metadata_id} \
			--node-data {input.branch_lengths} {input.traits} {input.nucleotide_mutations} {input.amino_acid_mutations} \
			--auspice-config {input.source_configuration} {input.types_configuration} \
			--title {params.title} \
			--maintainers {params.maintainers} \
			--build-url {params.build_url} \
			--description {input.description} \
			--colors {input.colors} \
			--lat-longs {input.coordinates} \
			--output {output} \
			--include-root-sequence-inline \
			--minify-json
		"""
