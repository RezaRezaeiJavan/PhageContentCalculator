##########General Information##########
# PhageContentCalculator Version 1.0
# 
# This script is designed to estimate the prophage content within bacterial genomes 
# 
# Author: Reza Rezaei Javan
# Copyright (C) 2019 University of Oxford
# E-mail: Reza.RezaeiJavan@ndm.ox.ac.uk

##########License Information##########
# This is free software; you can redistribute it under
# the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
#  This is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY.See the GNU Lesser General
#  Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this software package; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA

##########VERY IMPORTANT!!!!!##########
##########Paths to Roary and Prokka - MUST BE DEFINED BY USER##########
#This depends on how you installed the programme and what version you are using. For example, this could be something like "perl /home/reza/prokka-1.11/bin/prokka" or just "prokka"
PATH_TO_PROKKA = "perl /home/ubuntu/prokka/bin/prokka"
PATH_TO_ROARY = "roary"

##########Import Packages##########
import csv
import os.path
import re
import shutil
import glob
import subprocess
import sys
import os

##########Parameters##########
#defining input
Input_file = sys.argv[1]
#threshold to use for roary
roary_threshold = 70

##########Functions##########
#Defining moving tagged files into a mixed group
def copyreferencefiles(path,group):  
    # retrieve file list
    filelist=glob.glob(path)
    for single_file in filelist:
         # move file with full paths as shutil.copy() parameters
        shutil.copy(single_file,"tmp/Inputs/{}/".format(group))

#check to see if folder exsist, if not create it. If it exists, remove it and recreate it again
def checkifexistifnotcreateit(path):
    newpath = r'{}'.format(path) 
    if os.path.exists(newpath):
        shutil.rmtree(newpath)
    if not os.path.exists(newpath):
        os.makedirs(newpath)

#Defining the use of roary
#You can add -p to use multiple CPUs, for example "-p 8" will use 8 CPUs.        
def roary(identity):  
    os.system('{} -i {} -f tmp/roary_results tmp/tagged_files_mixed/*.gff 2>/dev/null'.format(PATH_TO_ROARY, identity))

#Defining Tagging files with group names
def tagfiles(path, groupname):  
  # Open a file for writing and create it if it doesn't exist
    for fname in glob.glob("{}".format(path)):
        string = open(fname).read()
        new_str = re.sub("ID=", "ID={}".format(groupname), string)
        open(fname, 'w').write(new_str)

#Defining moving tagged files into a mixed group
def movefilesmixfolder(path):  
    # retrieve file list
    filelist=glob.glob(path)
    for single_file in filelist:
         # move file with full paths as shutil.move() parameters
        shutil.move(single_file,"tmp/tagged_files_mixed/")

#Annotating the input file
def annotating_input():
	#Create a working folder
	checkifexistifnotcreateit("tmp/Annotating_the_input")
	#Move the genome of interest to a specific folder in order to annotate it
	shutil.copy(Input_file,"tmp/Annotating_the_input")
	#Shorten file names if they are too long (to avoid Prokka crashing; Prokka won't annotate files that are too long)    
	filelist=glob.glob("tmp/Annotating_the_input/*.fas")
	for single_file in filelist:
	    contig_name_with_extension = single_file.split("/")[2]
	    contig_name = contig_name_with_extension.split(".")[0]
	    if len(contig_name) > contig_len_threshold:
	    	shortened_single_file = "tmp/Annotating_the_input/{}.fas".format(contig_name[0:9])
	    	os.rename(single_file, shortened_single_file)
	#Annotate the file using prokka (there should only be one file, but the script actually runs for as many file as there is)
	filelist=glob.glob("tmp/Annotating_the_input/*.fas")    	
	for single_file in filelist:
	    contig_name_with_extension = single_file.split("/")[2]
	    contig_name = contig_name_with_extension.split(".")[0]
	    os.system("{} --outdir tmp/Annotating_the_input/annotated/{}.prokka --kingdom Bacteria --locustag {} {} --centre XXX --compliant --force 2>/dev/null".format(PATH_TO_PROKKA, contig_name,contig_name,single_file))
	#move the annotated gff files to GroupC and delete the working folder
	for name in glob.glob('tmp/Annotating_the_input/annotated/*/*.gff'):
	    shutil.move(name,"tmp/Inputs/GroupC/")
	# shutil.rmtree("tmp/Annotating_the_input/")
#Tag files
def tag_files_function():
	if os.path.exists("tmp/Inputs/GroupA"):
	    tagfiles("tmp/Inputs/GroupA/*.gff","GroupA")
	if os.path.exists("tmp/Inputs/GroupB"):
	    tagfiles("tmp/Inputs/GroupB/*.gff","GroupB")
	if os.path.exists("tmp/Inputs/GroupC"):
	    tagfiles("tmp/Inputs/GroupC/*.gff","GroupC")

#Veendiagram
def veendiagram():
	if (os.path.isfile('GroupA.CSV')is True) and (os.path.isfile('GroupB.CSV') is True):
	    A = open('GroupA.CSV')

	    CSV_A = csv.reader(A)

	    GroupA = set([])
	    for row in CSV_A:
	        GroupA.add(row[0])

	    A.close()

	    B = open('GroupB.CSV')

	    CSV_B = csv.reader(B)

	    GroupB = set([])
	    for row in CSV_B:
	        GroupB.add(row[0])

	    B.close()
	else: print('**********Job aborted - Could not extract ORFs from the full-length and satellite prophage databases**********')

	#Import 3rd column from CSV file for GroupC
	if os.path.isfile('GroupC.CSV')is True:
	    C = open('GroupC.CSV')
	    CSV_C = csv.reader(C)
	    GroupC = set([])
	    for row in CSV_C:
	        GroupC.add(row[0])
	    C.close()
	else: pass


	###CALCULATIONS
	if (os.path.isfile('GroupA.CSV') is True and os.path.isfile('GroupB.CSV') is True and os.path.isfile('GroupC.CSV') is True):

	    #First round of calculations
	    A_in_B = GroupA.intersection(GroupB)
	    A_in_C = GroupA.intersection(GroupC)
	    A_in_B_in_C = A_in_B.intersection(GroupC)
	    A_in_B_not_C = A_in_B.difference(GroupC)
	    A_in_C_not_B = A_in_C.difference(GroupB)
	    Not_A = A_in_B.union(A_in_C_not_B)
	    A_only = GroupA.difference(Not_A)


	    #Second round of calculations
	    B_in_A = GroupB.intersection(GroupA)
	    B_in_C = GroupB.intersection(GroupC)
	    B_in_A_in_C = B_in_A.intersection(GroupC)
	    B_in_A_not_C = B_in_A.difference(GroupC)
	    B_in_C_not_A = B_in_C.difference(GroupA)
	    Not_B = B_in_A.union(B_in_C_not_A)
	    B_only = GroupB.difference(Not_B)


	    #Third round of calculations
	    C_in_A = GroupC.intersection(GroupA)
	    C_in_B = GroupC.intersection(GroupB)
	    C_in_A_in_B = C_in_A.intersection(GroupB)
	    C_in_A_not_B = C_in_A.difference(GroupB)
	    C_in_B_not_A = B_in_C.difference(GroupA)
	    Not_C = C_in_A.union(C_in_B_not_A)
	    C_only = GroupC.difference(Not_C)


	    ##PRINTING
	    input_name = Input_file.split(".")[0]
	    Total_gene_num = float(len(GroupC))
	    FP_origin_num = float(len(C_in_A))
	    SP_origin_num = float(len(C_in_B))
	    Unique_gene_num = float(len(C_only))

	    print("{} --- Full-length Phage Origin: {:.1%} --- Satellite Phage Origin: {:.1%}".format(input_name, float(FP_origin_num/Total_gene_num),float(SP_origin_num/Total_gene_num)))
	    
	    with open("results.txt", "a") as text_file:
	    	text_file.write("{}\t{:.3%}\t{:.3%}\t{}\t{}\t{}\t{}\n".format(input_name, float(FP_origin_num/Total_gene_num), float(SP_origin_num/Total_gene_num), int(Total_gene_num), int(FP_origin_num), int(SP_origin_num), int(Unique_gene_num)))

	else:
		print "WARNING: could not do the calculations."

#delete nnnecessary files to tidy things up for the next round
def delete_unnecessary_files():
	shutil.rmtree("tmp/")
	os.remove("GroupA.CSV")
	os.remove("GroupB.CSV")
	os.remove("GroupC.CSV")

##########Tidying up files##########
###Tidying up files to ensure they are compatable with Prokka###
#shortening contigs above this number of characters to ensure compatibility with Prokka
contig_len_threshold = 10


##START ANALYSIS###
#check to see if folder exsist, if not create it
try:
	checkifexistifnotcreateit("tmp/Inputs/GroupA")
	checkifexistifnotcreateit("tmp/Inputs/GroupB")
	checkifexistifnotcreateit("tmp/Inputs/GroupC")

	#Copy full-length phage and satellite phage files into Inputs
	copyreferencefiles("db/FP/*.gff","GroupA")
	copyreferencefiles("db/SP/*.gff","GroupB")

	try:
		try:
			#Annotate the input
			annotating_input()
		except:
			print "ERROR: something went wrong with annotating the input"

		try:
			#tag files
			tag_files_function()
		except:
			print "ERROR: something went wrong with tagging the files"

		try:
			##Mixing all these files into one folder
			#check to see if folder exsist, if not create it
			checkifexistifnotcreateit("tmp/tagged_files_mixed")
			#move them into the mix folder
			movefilesmixfolder("tmp/Inputs/GroupA/*.gff")
			movefilesmixfolder("tmp/Inputs/GroupB/*.gff")
			movefilesmixfolder("tmp/Inputs/GroupC/*.gff")
		except:
			print "ERROR: something went wrong with mixing the files to do Roary on"

		try:
			###Running Roary
			roary(roary_threshold)
		except:
			print "ERROR: something went wronng with the Roary analysis"   

		try:
			#greping
			os.system("grep -e \"GroupA\" tmp/roary_results/gene_presence_absence.csv > GroupA.CSV")
			os.system("grep -e \"GroupB\" tmp/roary_results/gene_presence_absence.csv > GroupB.CSV")
			os.system("grep -e \"GroupC\" tmp/roary_results/gene_presence_absence.csv > GroupC.CSV")
		except:
			print "ERROR: something went wrong with grep"

		try:
			veendiagram()
		except:
			print "ERROR: something went wrong with the (veen diagram) calculations"
	except:
		print "Unknown ERROR... Something went wrong..."
	try:
		delete_unnecessary_files()
	except:
		print "ERROR: something went wrong with deleting unnecessary files"
except:
	print "Analysis failed completely... Could not even create the folders or remove unnecessary files "

