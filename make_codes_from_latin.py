#!/usr/bin/env python3
#
import argparse


__author__ = "Ekaterina Osipova, 2022."


## Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, help='input file with latin names')
args = parser.parse_args()

## Read file Apus apus tree -> convert to apuApuTre
with open(args.input, 'r') as inf:
	for line in inf.readlines():
		components = line.split()
		if len(components) == 1:
			components = line.split('_')
			if len(components) == 1:
				print('could not convertline: {}'.format(line))
		
		c = components[0]
		code = c[:3].lower()
		for c in components[1:]:
			code += c[0].upper() + c[1:3]
		print(code)
