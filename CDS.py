#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 18:32:37 2020

@author: gin
"""
import numpy as np

# Script with functions to read/ manipulate CDS data

def get_type(letter):
    '''
    Returns function to convert string to relevant numerical type 
    (CDS chooses from 'I', 'F', 'A' = int, float, str)
    Example: 
        >>> get_type('F')('3')
        >>> 3.0
    '''
    if letter == 'I':
        return int
    if letter == 'F':
        return np.float64
    if letter == 'A':
        return str

    
def read_table(readmefile, tablefile):
    
    # Dictionary that holds column name+byte info
    column_dict = {}
    
    
    fop = open(readmefile)
    start_read = False
    end_read = False
    jj = 0 # will count two lines after the start point
    for line in fop.readlines():
        
        # End statement: if you're done reading lines, exit this loop
        if end_read ==True:
            break
        
        # Statement that prepares for exiting the loop
        # The column info stars after the first line satisfying this check
        if line[0] == '-' and start_read== True:
            jj+=1
            
            # The column info ends at the second line in the file that satisfied 
            # the above check. 
            if jj > 1:
                end_read = True
                continue
            else:
                continue
            
        # Read the column names and byte info
        if start_read == True:
            
            # First check that the line you'll read contains info and is not just
            # a continuation of the comment from the line above
            if line[4]==' ':
                continue
            
            # Create a dictionary entry named after the column that holds: 
            # first byte of data value, last byte, data type, empty list to be filled with data afterwards
            
            # find start of each column in the Readme file explanation
            sl = line.split()
            sl2 = line.split('-') 
            
            # first non-space group of characters 
            start_idx = line.find(sl[0])
            # second non-space group of characters
            end_idx = line.find(sl[1])
            end_idx2 = line.find(sl2[1]) # second case if no gap between start and end byte values
            # third
            format_idx = line.find(sl[2])
            format_idx2 = line.find(sl[1])
            # fifth
            label_idx  = line.find(sl[4])
            label_idx2 = line.find(sl[3])
            # final
            explans_idx = line.find(sl[5])
            explans_idx2 = line.find(sl[4])
            
            #column_dict[line[21:28].strip()] = [int(line[0:4]), int(line[5:8]), get_type(line[9]), []]
            try:
                start_byte = int(line[start_idx:end_idx].strip(' ').strip('-')) # remove trailing '-'
                end_byte = int(line[end_idx:format_idx].strip())
                label = line[label_idx:explans_idx].strip()
                
                column_dict[label] = [start_byte, end_byte, get_type(line[format_idx]), []]
                
            except ValueError:
                start_byte = int(line[start_idx:start_idx+3].strip('-'))
                end_byte = int(line[end_idx2:format_idx2].strip())
                label = line[label_idx2:explans_idx2].strip()
            
                column_dict[label] = [start_byte, end_byte, get_type(line[format_idx2]), []]
            
        # Here is where we find which line is the marker to start reading in the column info
        sl = line.split()
        # These words appear two lines before the info on table columns
        if 'Bytes'  in sl and 'Format' in sl:
            # Set marker to True
            start_read = True
        
            
    fop.close()

    
    # Now read the table file
    fop = open(tablefile)
    for nline in fop.readlines():
        # skip empty lines
        if nline == '\n':
            continue
        line = nline.strip('\n')
        # Loop over columns (specified in column_dict.keys())
        for key in column_dict.keys():
            # get the start and end byte of the value we need to read
            cstart = column_dict[key][0]-1 # python numbering
            cend = column_dict[key][1]
            # get the data type
            datatype = column_dict[key][2]
            
            # append data from this row to the given column
            column_dict[key][3].append(datatype(line[cstart:cend]))

    
    fop.close()        
    # Convert list of data in each column to array --  quick fix! We'll have to
    # move to a better format, e.g. astropy tables or VOTable
    for key in column_dict.keys():
        column_dict[key][3] = np.array(column_dict[key][3])
    
    return column_dict
        
        
        
# Run example
coldict = read_table('example_CDS_README.txt','example_CDS_table.txt')       
            
            
        
