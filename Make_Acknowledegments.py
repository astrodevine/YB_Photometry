# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 11:53:47 2024

@author: colemaya
"""

import pandas as pd
import os
from operator import itemgetter
from datetime import date
#type 'pip install better-profanity' in terminal for this code to work
from better_profanity import profanity

##Functions##

def sort_by_last_name(name_list):
    # sort the list of names by the last name using itemgetter
    name_list = [x.split() for x in name_list]
    sorted_list = [" ".join(i)  for i in sorted(name_list, key=itemgetter(-1))]
    return sorted_list

def format_acknowledgements(name_list):
    # function to clean, format, sort, and filter the list
    cleanednamelist = [x for x in name_list if str(x) != 'nan']
    cleanednamelist = list(set(cleanednamelist))
    #capitalize first and last names and remove extra spaces
    formattednames= []
    for name in cleanednamelist:
        newname = name.strip()
        elem = newname.split(' ')
        elem[0] = elem[0][0].upper() + elem[0][1:]
        elem[-1] = elem[-1][0].upper() + elem[-1][1:]
        #exclude submissions with only a first name
        if len(elem) > 1:
            formattednames.append(' '.join(elem))
    #remove Dr. Devine's fake example entry
    formattednames.remove('Ima Q. Student')
    censorednames = remove_profanity(formattednames)
    sortednames = sort_by_last_name(censorednames)
    return sortednames


def remove_profanity(name_list):
    # identifies potentially inappropriate names and writes them to a list for manual review. user presses enter to finish
    flag_list= []
    for name in name_list:
        if profanity.contains_profanity(name) or name.replace(' ','').replace('.','').replace('-','').isalpha() == False:
            flag_list.append(name)
    removeindex = 0
    while removeindex in range(len(flag_list)):
        print('Potentially Problematic Names:')
        for name in flag_list:
            print(str(flag_list.index(name) + 1) + ': ' + name)
        removeinput = input('Enter the number corresponding to the name you would like to remove, or press enter when done:')
        if removeinput.isnumeric():
            removeindex = int(removeinput) - 1
            if removeindex in range(len(flag_list)):
                name_list.remove(flag_list[removeindex])
        else:
            removeindex = -1
    return name_list

##Code##

#retrieve list of names from file
path = '.'
ack_filename = os.path.join(path, 'YB_Acknowledgements.csv') #edit this line to match your csv file name
ack_file = pd.read_csv(ack_filename)
currentheaders = ack_file.columns.tolist()
names = ack_file[currentheaders[1]]
namelist= names.tolist()

nicelist = ['Contributors'] + format_acknowledgements(namelist)

#write to csv file labeled with date created
with open('YB_acklist_' + str(date.today()) + '.csv', 'w+') as f: 
    # write elements of list
    for items in nicelist:
        f.write('%s\n' %items)
 
# close the file
f.close()

    
    
    
    
    
