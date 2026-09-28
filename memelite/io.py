# io.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import re
import numpy


def read_meme(filename, n_motifs=None):
	"""Read a MEME file and return a dictionary of PWMs.

	This method takes in the filename of a MEME-formatted file to read in
	and returns a dictionary of the PWMs where the keys are the metadata
	line and the values are the PWMs. Each key is the rest of the motif's
	MOTIF line, such as 'MA0004.1 Arnt', with its fields joined by one space.


	Parameters
	----------
	filename: str
		The filename of the MEME-formatted file to read in


	Returns
	-------
	motifs: dict
		A dictionary of the motifs in the MEME file.
	"""

	motifs = {}

	# The limit is checked after each motif is added, which never stops at 0.
	if n_motifs == 0:
		return motifs

	with open(filename, "r") as infile:
		motif, width, i = None, None, 0

		for line in infile:
			if motif is None:
				if line[:5] == 'MOTIF':
					# The fields may be separated by tabs or several spaces.
					motif = ' '.join(line[5:].split())
				else:
					continue

			elif width is None:
				if line[:6] == 'letter':
					width = int(re.search(r'\bw=\s*(\d+)', line).group(1))
					pwm = numpy.zeros((width, 4))

			else:
				pwm[i] = list(map(float, line.strip("\r\n").split()))
				i += 1

				# Stored as soon as the last row is read, rather than on the
				# line after it, which may be the next MOTIF line or absent.
				if i == width:
					motifs[motif] = pwm.T
					motif, width, i = None, None, 0

					if n_motifs is not None and len(motifs) == n_motifs:
						break

	return motifs


def write_meme(filename, motifs):
	"""Write a MEME file.

	This method takes in a filename and either a list or dictionary of motifs and
	writes them to disk in a MEME-formatted file.


	Parameters
	----------
	filename: str
		The name of the MEME-formatted file to save.

	motifs: list or dict
		The set of motifs to save. If a list, the name of each motif will be its
		numerical ordering in the list. If a dictionary, the name will be the key
		in the dictionary.
	"""

	with open(filename, "w") as outfile:
		outfile.write("MEME version 4\n\n")
		outfile.write("ALPHABET= ACGT\n\n")
		outfile.write("strands: + -\n\n")
		outfile.write("Background letter frequencies\n")
		outfile.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")

		if isinstance(motifs, dict):
			motif_pwms = list(motifs.values())
			motif_names = list(motifs.keys())
		else:
			motif_pwms = motifs
			motif_names = [str(i) for i in range(len(motifs))]

		for name, pwm in zip(motif_names, motif_pwms):
			outfile.write("MOTIF {}\n".format(name))
			outfile.write("letter-probability matrix: alength= {} w= {} nsites= 1 E= 0\n".format(*pwm.shape))

			for col in pwm.T:
				outfile.write("{} {} {} {}\n".format(*col))

			outfile.write("URL BLANK\n\n")
		