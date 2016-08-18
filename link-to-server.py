import os,sys
species = sys.argv[1]
def link_data_to_server():
	"""
	Sending json files and alignment to server for data visualization
	"""
	os.system('rsync -avz ./data/%s/meta-dict-%s.js root@mpmtest.eb.local:/root/server-side/public/javascripts/ '%(species,species))
	os.system('rsync -avz ./data/%s/vis/ root@mpmtest.eb.local:/root/server-side/public/dataset/%s/ '%(species,species))
link_data_to_server()