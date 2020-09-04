#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metagenomic Mock Community Generator (MMCG) v0.1.0
Uses wgsim, nanosim and seqtk to generate mock communities
Comma-delimiter separated file used to specify the proportions of the metagenome
"""
   
import subprocess
import os
import datetime
import time
import argparse
import numpy

"""
Function for obtaining the current time
"""
def curr_time():
   ts = time.time()
   st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
   return(st)
   
   

"""



Function for obtaining directories at command line
"""
def dir_path(string):
    real_path=os.path.realpath(string)
    if not os.path.isdir(string):
        print("Creating output directory: {}".format(real_path))
        os.makedirs(real_path)
    return real_path
"""
Command-line arguments
"""

parser = argparse.ArgumentParser(description='Generate mock metagenomic read data.', formatter_class=argparse.RawTextHelpFormatter)
required_arg = parser.add_argument_group('required named arguments')
required_arg.add_argument('-i', metavar ='', type=str, nargs=1, required=True,
                   help='Input abundance file. Comma-delimited file containing file directory, total read number, file names and abundances')
required_arg.add_argument('-o', metavar ='', type=str, nargs=1, required=True,
                   help='''Output name. Read output = <output>_1.fastq and <output>_2.fastq. Abundance output = <output>_wgsim_reads.csv
Warning. Do not supply path name as all output will be supplied to a default or given (see -o_dir) directory''')
required_arg.add_argument('-s', metavar ='', type=str, nargs=1, choices = ['miseq', 'ont'], required=True,
                   help='''Sequencer choice. Options (must select one from key): 
miseq = Illumina MiSeq, paired-end, 300bp, error rate = 0.2%%
ont = Oxford Nanopore, long-read, error rate = 5-20%%''')

optional_arg = parser.add_argument_group('optional name arguments')
optional_arg.add_argument('-o_dir', type=dir_path, nargs=1,
                   help='''Output directory. 
Default: Generate new folder (mmcg_out) in current working directory''')

args = parser.parse_args()

# Obtain current directories: terminal directory, relative script directory
terminal_cwd = os.getcwd()
output_directory = terminal_cwd
script_cwd = os.path.dirname(os.path.realpath(__file__))

# Print the start of the Python script
print("Start time: " + curr_time())

# Create output directory
# Check if optional argument has been supplied and create 
if not args.o_dir:
    args.o_dir = os.path.realpath("mmcg_out")
    if not os.path.isdir(args.o_dir):
        print("Creating output directory: {}".format(args.o_dir))
        os.makedirs(args.o_dir)
        
# Save command-line names
input_file_name = os.path.realpath(args.i[0])
seq_choice = args.s[0]
output_name = args.o[0]
output_dir = args.o_dir
# Checks if the args.o_dir is in list format and converts it to string
if isinstance(output_dir, list):
    output_dir = output_dir[0]

# Check that input file exists
if not os.path.isfile(input_file_name):
    print("Input file/path does not exist")
    print("End time: " + curr_time())
    exit()
    
# Add slash for output directory if none is given
if not output_dir[-1] == "/":
    output_dir = output_dir + "/"
    

file_delim = ","
# Open and read file containing file names, abundances, total number of reads and directory
input_file = open(input_file_name)
input_lines = input_file.readlines()
input_file.close()



# Save file information
file_info = {}
# Extract directory and total number of reads
file_dir = input_lines[0].split(file_delim)[0]
if not file_dir[-1] == "/":
    file_dir = file_dir + "/"
total_number_reads = int(input_lines[1].split(file_delim)[0])

# TEST
print("\nTEST")
print("Total number of reads: {}".format(total_number_reads))



# Remove these lines, so that only file name and abundances are present
del input_lines[:2]

# TEST
print("\nTEST")
print("Usable lines: {}".format(input_lines))

# Parameter for sequencing simulation
# This can be expanded upon to run different sequencers
parameter_dict = {'miseq': ['300', '300', '0.002'],
                  'ont' : 0}

seq_parameters = parameter_dict[seq_choice]

#TEST
print("\nTEST")
print("Sequencing parameters: {}".format(seq_parameters))


# Open summary file to write read simulation results
read_overview_filename = output_dir + output_name + "_overview.csv"
read_overview_file=open(read_overview_filename, "w")

#line =input_lines[0]
for line in input_lines:
    print("Run start time: " + curr_time())
    line_split = line.rstrip("\n").split(file_delim)
    
    # Extract the ID for each file to add to simulated reads
    acc_ID = line_split[0].split("_")[0] + "_" + line_split[0].split("_")[1]
    
    
    # Combine the directory with filename
    genome_file = file_dir + line_split[0]
    
    
    # Increase the read number which can be subsampled later to ensure the correct number of reads are generated
    sim_read_number = int(total_number_reads * (float(line_split[1]) + 0.1))
    seqkt_read_number = int(total_number_reads * (float(line_split[1])))
    
    
    # If the sequencer selected is a short-read sequencer it will have parameters stored as an array
    # This can be used to select which sequencer is used
    if isinstance(seq_parameters, list):
        # Pass parameters into "wgsim.sh" to perform short-read sequencing
        # Create an input list
        os.chdir(script_cwd)
        wgsim_input_raw=['bash',"wgsim.sh", seq_parameters[0], 
                         seq_parameters[1], seq_parameters[2], 
                            sim_read_number, seqkt_read_number,
                            genome_file, output_name, acc_ID, output_dir]
        
        # Convert all list elements to strings
        wgsim_input=[str(elem) for elem in wgsim_input_raw]
        
        # Run wgsim
        wgsim_process=subprocess.run(wgsim_input, stdout=subprocess.PIPE)
        wgsim_out=wgsim_process.stdout.decode().splitlines()
        print(wgsim_out)
        
        # Store and write abundances
        wgsim_abundance_1 = str(int(wgsim_out[0])/total_number_reads)
        wgsim_abundance_2 = str(int(wgsim_out[1])/total_number_reads)
        
        seqkt_abundance_1 = str(int(wgsim_out[2])/total_number_reads)
        seqkt_abundance_2 = str(int(wgsim_out[3])/total_number_reads)
        
        wgsim_read_output = [acc_ID] + wgsim_out[0:2] + \
        [wgsim_abundance_1, wgsim_abundance_2] + wgsim_out[2:] + [seqkt_abundance_1,
        seqkt_abundance_2]    
        read_overview = ",".join(wgsim_read_output)
        read_overview_file.write(read_overview + "\n") 
    elif isinstance(seq_parameters, int):
        # Generate nanopore reads
        if seq_parameters == 0:
            # Change directory to the nanosim model
            os.chdir(script_cwd+"/nanosim_model")
    
            
            # Pass parameters into "nanosim.sh" to perform short-read sequencing
            # Create an input list
            nanosim_input_raw=['bash',"nanosim.sh", 
                                sim_read_number, seqkt_read_number,
                                genome_file, output_name, acc_ID, output_dir]
            
            # Convert all list elements to strings
            nanosim_input=[str(elem) for elem in nanosim_input_raw]
            
            # Run nanosim
            nanosim_process=subprocess.run(nanosim_input, stdout=subprocess.PIPE)
            nanosim_out=nanosim_process.stdout.decode().splitlines()
            print(nanosim_out)
            
            # Store and write abundances
            nanosim_abundance_1 = str(int(nanosim_out[0])/total_number_reads)
             
            seqkt_abundance_1 = str(int(nanosim_out[1])/total_number_reads)
            
            nanosim_read_output = [acc_ID] + [nanosim_out[0]] + \
            [nanosim_abundance_1] + [nanosim_out[1]] + [seqkt_abundance_1]    
            read_overview = ",".join(nanosim_read_output)
            read_overview_file.write(read_overview + "\n")
            

# Change directory to script working directory
os.chdir(script_cwd)

