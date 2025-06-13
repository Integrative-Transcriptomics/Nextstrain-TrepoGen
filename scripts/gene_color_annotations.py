#!/usr/bin/env python3
"""Color annotations in a JSON file produced by augur translate.

Description:
	This script updates the color annotations in a JSON file produced by augur translate,
	using a fixed color map for genome annotations (based on the managua color map).
	It reads the input JSON file, updates the "color" field of each annotation whose name
	matches a key in the predefined color map, and writes the modified data back to the file.

Usage:
	python gene_color_annotations.py -f <input_json_file>
"""
import json
import argparse

"""Fixed color map for genome annotations in Auspice JSON files; based on the managua color map (https://www.fabiocrameri.ch/colourmaps/)."""
color_map = {
	"Signal-peptide": "#74c4ea",
	"Hatch": "#68a3d5",
	"Intracellular-loop": "#5d85c0",
	"Transmembrane-helix": "#893e3a",
	"Beta-strand": "#4d4c88",
	"Extracellular-loop": "#e9ac5b"
}

def parse_args():
	"""Parses command-line arguments for the gene color annotation script.

	Returns:
		argparse.Namespace: Parsed command-line arguments with the following attributes:
			file (str): Path to the input file (required).
	"""
	parser = argparse.ArgumentParser(
		description="",
	)
	parser.add_argument("-f", "--file", type=str, required=True, help="")
	return parser.parse_args()

def main():
	"""Set colors for annotations in a JSON file produced by augur translate.

	Description:
		Reads a JSON file produced by augur translate, updates the color of each annotation
		based on a predefined color map, and writes the modified data to the same file.
		For each annotation in the input, if its name starts with a key from the global
		color_map, its "color" field is updated accordingly.

	Arguments:
		None. Arguments are parsed internally via parse_args().

	Returns:
		None.

	Side Effects:
		The input JSON file is modified in place, with the "color" field of annotations updated.

	Raises:
		FileNotFoundError: If the input file does not exist.
		json.JSONDecodeError: If the input file is not valid JSON.
		KeyError: If expected keys are missing in the JSON structure.
	"""
	args = parse_args()

	global color_map

	data = {}
	with open(args.file, 'r') as i:
		data = json.load(i)

	if "annotations" in data:
		for annotation in data["annotations"].keys():
			for color_key, color_value in color_map.items():
				if annotation.startswith(color_key):
					data["annotations"][annotation]["color"] = color_value
					break

	with open(args.file, 'w+') as o:
		json.dump(data, o, indent=2)

if __name__ == '__main__':
	main()
