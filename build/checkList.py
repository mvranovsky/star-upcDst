#!/usr/bin/python



input2 = "/star/u/mvranovsk/star-upcDst/build/input.list"


with open(input2, 'r') as source2_file, open("inputcheck.list", 'w') as dest_file:
	for line in source2_file:
		reduced = line.replace("/star/u/adamczyk/pwd/picoUPC/test_upcdst2/","").split("/")
		#print reduced
		line2 = "/star/data01/pwg_tasks/upc03/pp17/UPC_P20ic/" + reduced[1]
		dest_file.write(line2)



