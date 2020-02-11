#!/usr/local/bin/python3
#
# Program will take the output from getting the COG categories and
# positions (not accounting for multiple COG cats) and will change it
# so that each gene will have one poition and either a 0 or a 1 for
# each of the COG categories. this will be used as input for R to
# calculate the logistic regression.
#
# run as "change_COG_to_0_1.py *_COG_category_and_positions.txt" 
#
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
import pandas as pd #importing pandas
# FOR DEBUGGING ONLY pdb.set_trace()
#
#below will open the file using pandas which apparently looks at stuff
#like a matrix

data = pd.read_csv(sys.argv[1], sep="\t", header = None)
data.columns = ["COG", "start", "end"]

#finding midpoint of each gene and adding it to new column in
#dataframe
data['midpoint'] = (data['start'] + data['end']) / 2

#rounding the midpoint column so it does not have any decimals
data.midpoint = data.midpoint.round()

#creat vector with all the COG categories
COGs = ['J', 'A', 'K', 'L', 'B', 'D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O', 'X', 'C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q', 'R', 'S']

#create empty vectors for each column
COG = []
midpoint = []
value = []

#go through each row in the dataframe
for index, row in data.iterrows():
#	print (row.COG.str.contains("[A-Z]{2,}"))
#if there is more than one COG per gene
	if (re.search("[A-Z]{2,}",row.COG)):
		mult_cogs = row.COG
		mult_cogs_list = list(mult_cogs)
		num_cogs = len(mult_cogs_list)
		for cog in COGs:
			for mult_cog in mult_cogs_list:
				if (mult_cog == cog):
					present = 1
					break
				else:
					present = 0
#if the mult cog is the same as the cog
			if (present == 1):
				COG.extend(cog)	
				mid_pt = row.midpoint
				midpoint.extend([mid_pt])	
				value.extend([1])	
			else:
				COG.extend(cog)	
				mid_pt = row.midpoint
				midpoint.extend([mid_pt])	
				value.extend([0])	
			present = None
#if there is one COG per gene
	else:
		for cog in COGs:
			if (row.COG == cog):
#adding values of 0 or 1 for each COG cat for each gene. all have the
#same midpoint. making each column of the dataframe separatly					
				COG.extend(cog)	
				mid_pt = row.midpoint
				midpoint.extend([mid_pt])	
				value.extend([1])	
			else:
				COG.extend(cog)	
				mid_pt = row.midpoint
				midpoint.extend([mid_pt])	
				value.extend([0])	

#combine the above columns into a data frame with columns named
new_cog_dat = pd.DataFrame({'COG': COG, 'midpoint': midpoint, 'value':
value})


#save data frame file to csv
new_cog_dat.to_csv('COG_binary_data.csv')

