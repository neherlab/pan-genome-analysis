#!/usr/bin/env python
import os,sys
import argparse

parser = argparse.ArgumentParser(description=\
    'Sending panX analysis results (clusters, alignment, tree, etc) to the internal/external server for pan-genome visualization and exploration',
    usage='./%(prog)s'+' -s E_coli -f /usr/pan-genome-visualization (-h: help)')
parser.add_argument('-s', '--species_folder', type = str, required=True,
    help='the folder where panX analysis pipeline runs, e.g.: E_coli (if the folder path is like this: /yourPath/pangenome_analysis/data/E_coli)', metavar='')
parser.add_argument('-f', '--visualization_folder', type = str, required=True,
    help='the absolute path for the visualization folder on the server, e.g.: /yourServer/pan-genome-visualization ', metavar='')

params  = parser.parse_args()
species = params.species_folder
vis_folder = params.visualization_folder

os.system('rsync -az ./data/%s/vis/ %s/public/dataset/%s/ '%(species,vis_folder,species))
# for external sever:
# 'rsync -az ./data/E_coli_1000/vis/ admin@xyz.abc.local:/usr/pan-genome-visualization/public/dataset/E_coli_1000/