import os, sys
import pandas as pd
import numpy as np
from collections import OrderedDict


def Parse_GFF(GFF, path_GFF, path_VR):
	os.chdir(path_GFF)

	out_f = GFF.replace('.gtf', '_biotypes.txt')
	##Scaffold\tGeneID\tTranscriptID\t\ProteinID\tProductID\tProduct\n


	list_GeneID=[]
	gff = open(GFF, 'r')
	lines = gff.readlines()


	scf=''
	GeneID=''
	BioType=''
	out_df = pd.DataFrame(columns=['Scaffold', 'GeneID', 'BioType'])
	for l in lines:
		if '#' not in l:
			l_split = l.split('\t')
			if l_split[2] == 'gene':
				scf = l_split[0]
				infos = l_split[8]

				infos_split = infos.replace('\n', '').split('; ')
				for i in infos_split:
					i_split = i.split(' ')
					if i_split[0] == 'gene':
						GeneID=i_split[1].replace('"', '')

					if i_split[0] == 'gene_biotype':
						BioType=i_split[1].replace('"', '')

				tmp_df = {'Scaffold': scf, 'GeneID': GeneID, 'BioType': BioType}
				out_df = out_df.append(tmp_df, ignore_index=True)

				GeneID=''
				BioType=''



	os.chdir(path_VR)
	out_df.to_csv(out_f, columns=out_df.columns, sep='\t', index=False)


GFF = sys.argv[1]
path_GFF = sys.argv[2]
path_VR = sys.argv[3]
Parse_GFF(GFF, path_GFF, path_VR)
