import glob

files = glob.glob('*.fas') # I am assuming that you have fasta files with extension .fas. If not then change *.fas to your alignment file extension, such as *.nex, *.fasta or *.phy

with open('RYinfo.txt', 'w') as fp:
    for file in files:
        fp.write('%s, 3\n'%file)
