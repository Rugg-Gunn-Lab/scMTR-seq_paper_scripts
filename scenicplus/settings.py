import io
import numpy as np
import anndata as anndata
import scipy as s
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import csv
import os
from dfply import *
# from IPython.display import display
from re import search
from pathlib import Path
import argparse


#########
## I/O ##
#########

io = {}

host = os.uname()[1]
print(host)
if search("capstone|pebble|compute", host):
	io["basedir"] = "/bi/home/lij/data/projects/scCUT_Tag_rna/E45"
	io["gene_metadata"] = "/Users/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
elif search("ebi", host):
	io["basedir"] = "/hps/nobackup2/research/stegle/users/ricard/gastrulation_multiome_10x"
	io["gene_metadata"] = "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
elif search("Workstation", host):
    io["basedir"]="/home/lijingyu/gastrulation/data/gastrulation_multiome_10x"                   
    io["gene_metadata"]='/home/lijingyu/gastrulation/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt'

else:
	print("Computer not recognised"); exit()

io["metadata"] = io["basedir"] + "/sample_metadata.txt.gz"
io["anndata"] = io["basedir"] + "/processed/rna/anndata.h5ad"
io["anndata_scvelo"] = io["basedir"] + "/processed/rna/velocyto/anndata_scvelo.h5ad"

#############
## Options ##
#############

opts = {}

# opts["samples"] = [
# 	"E7.5_rep1",
# 	"E7.5_rep2",
# 	"E8.0_rep1",
# 	"E8.0_rep2",
# 	"E8.5_rep1",
# 	"E8.5_rep2"
# ]

# opts["rename_samples"] = {
#     "E8.5_rep1" : "multiome1",
#     "E8.5_rep2" : "multiome2",
#     "E7.5_rep1" : "rep1_L001_multiome",
#     "E7.5_rep2" : "rep2_L002_multiome"
# }


opts["stages"] = [
	"D0",
	"D1",
	"D2",
	'D3'
]
opts["samples"] = opts["stages"]
# opts["celltypes"] = [
# 	"Epiblast",
# 	"Primitive_Streak",
# 	"Caudal_epiblast",
# 	"PGC",
# 	"Anterior_Primitive_Streak",
# 	"Notochord",
# 	"Def._endoderm",
# 	"Gut",
# 	"Nascent_mesoderm",
# 	"Mixed_mesoderm",
# 	"Intermediate_mesoderm",
# 	"Caudal_Mesoderm",
# 	"Paraxial_mesoderm",
# 	"Somitic_mesoderm",
# 	"Pharyngeal_mesoderm",
# 	"Cardiomyocytes",
# 	"Allantois",
# 	"ExE_mesoderm",
# 	"Mesenchyme",
# 	"Haematoendothelial_progenitors",
# 	"Endothelium",
# 	"Blood_progenitors_1",
# 	"Blood_progenitors_2",
# 	"Erythroid1",
# 	"Erythroid2",
# 	"Erythroid3",
# 	"NMP",
# 	"Rostral_neurectoderm",
# 	"Caudal_neurectoderm",
# 	"Neural_crest",
# 	"Forebrain_Midbrain_Hindbrain",
# 	"Spinal_cord",
# 	"Surface_ectoderm",
# 	"Visceral_endoderm",
# 	"ExE_endoderm",
# 	"ExE_ectoderm",
# 	"Parietal_endoderm"
# ]

# opts["celltype_colors"] = {
# 	"Epiblast" : "#635547",
# 	"Primitive_Streak" : "#DABE99",
# 	"Caudal_epiblast" : "#9e6762",
# 	"PGC" : "#FACB12",
# 	"Anterior_Primitive_Streak" : "#c19f70",
# 	"Notochord" : "#0F4A9C",
# 	"Def._endoderm" : "#F397C0",
# 	"Gut" : "#EF5A9D",
# 	"Nascent_mesoderm" : "#C594BF",
# 	"Mixed_mesoderm" : "#DFCDE4",
# 	"Intermediate_mesoderm" : "#139992",
# 	"Caudal_Mesoderm" : "#3F84AA",
# 	"Paraxial_mesoderm" : "#8DB5CE",
# 	"Somitic_mesoderm" : "#005579",
# 	"Pharyngeal_mesoderm" : "#C9EBFB",
# 	"Cardiomyocytes" : "#B51D8D",
# 	"Allantois" : "#532C8A",
# 	"ExE_mesoderm" : "#8870ad",
# 	"Mesenchyme" : "#cc7818",
# 	"Haematoendothelial_progenitors" : "#FBBE92",
# 	"Endothelium" : "#ff891c",
# 	"Blood_progenitors" : "#c9a997",
# 	"Blood_progenitors_1" : "#f9decf",
# 	"Blood_progenitors_2" : "#c9a997",
# 	"Erythroid" : "#EF4E22",
# 	"Erythroid1" : "#C72228",
# 	"Erythroid2" : "#f79083",
# 	"Erythroid3" : "#EF4E22",
# 	"NMP" : "#8EC792",
# 	"Neurectoderm" : "#65A83E",
# 	"Rostral_neurectoderm" : "#65A83E",
# 	"Caudal_neurectoderm" : "#354E23",
# 	"Neural_crest" : "#C3C388",
# 	"Forebrain_Midbrain_Hindbrain" : "#647a4f",
# 	"Spinal_cord" : "#CDE088",
# 	"Surface_ectoderm" : "#f7f79e",
# 	"Visceral_endoderm" : "#F6BFCB",
# 	"ExE_endoderm" : "#7F6874",
# 	"ExE_ectoderm" : "#989898",
# 	"Parietal_endoderm" : "#1A1A1A"
# }

# opts["stages_colors"] = {
# 	'E6.5':"#D53E4F",
# 	'E6.75':"#F46D43",
# 	'E7.0':"#FDAE61",
# 	'E7.5':"#FFFFBF",
# 	'E7.25':"#FEE08B",
# 	'E7.75':"#E6F598",
# 	'E8.0':"#ABDDA4",
# 	'E8.5':"#3288BD",
# 	'E8.25':"#66C2A5",
# 	'mixed_gastrulation': "#A9A9A9"  
# }
