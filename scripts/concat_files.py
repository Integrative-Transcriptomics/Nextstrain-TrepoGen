#!/usr/bin/env python3
"""Merges multiple files into a single output file.

Description:
	If all input and output files have a `.json` extension, it merges the 
	JSON objects from each input file into a single dictionary and writes
	the result as JSON to the output file. The merge is flat, meaning only
	the first-level keys are considered and values are not replaced if a
	key already exists. However, if a key's value is a list, it extends
	that list with the values from subsequent files.

	For other file types, the script reads the contents of each input file in
	the order provided, concatenates them, and writes the combined content to
	the specified output file.

Usage:
	python concat_files.py <input1.ext> <input2.ext> ... -o <merged_output.ext>

Note:
	No further validation is executed on the file's contents.
"""
import argparse
import json

def parse_args():
	"""Parses command-line arguments for concatenating files.

	Description:
		Parses command-line arguments to obtain a list of input file paths and
		one output file path. No validation is performed on the file contents or
		their formats.

	Arguments:
		-o, --output (required): Path to the output merged file.
		inputs (positional, required): One or more paths to input files to merge.

	Returns:
		argparse.Namespace: Parsed command-line arguments with the following attributes:
		- inputs (list of str): List of input file paths.
		- output (str): Path to the output file.
	"""
	parser = argparse.ArgumentParser(
		description="Merges multiple files into a single output file."
	)
	parser.add_argument('inputs', nargs='+', help='Paths to files to merge.')
	parser.add_argument('-o', '--output', required=True, help='Path to the output merged file.')
	return parser.parse_args()

def main():
	"""Concatenates the contents of multiple input files into a single output file.

	Description:
		This function parses command-line arguments to obtain a list of input file paths
		and one output file path.

		1. If all input files and the output file have a `.json` extension, it merges the 
		JSON objects loaded from each input file into a single object and writes the result
		as one object to the JSON output file. This is done in a flat way, i.e. only the 
		keys of the first level of the input objects are considered - whereby array contents
		are extended - and values are not replaced.

		2. For other file types, it reads the contents of each input file, concatenates them
		in order, and writes the combined content to the specified output file.

	Arguments:
		None. Arguments are parsed internally via `parse_args()`.

	Returns:
		None.

	Side Effects:
		Writes the merged content to the output file specified in the CLI arguments.
	"""
	args = parse_args()
	if all( [ input.endswith( ".json" ) for input in args.inputs ] ) and args.output.endswith( ".json" ) :
		content = { }
		for input_file_path in args.inputs :
			with open( input_file_path, "r" ) as input_file :
				json_content = json.load( input_file )
				for key, value in json_content.items( ) :
					if type( value ) is list and key in content and type( content[ key ] ) is list :
						content[ key ].extend( value )
					else :
						content.setdefault( key, value )
		with open( args.output, "w+" ) as output_file :
			output_file.write( json.dumps( content, indent=2 ) )
	else :
		content = []
		for input_file_path in args.inputs :
			with open( input_file_path, "r" ) as input_file :
				content.extend( input_file.readlines( ) )
		with open( args.output, "w+" ) as output_file :
			output_file.write( "".join( content ) )

if __name__ == '__main__':
	main()
