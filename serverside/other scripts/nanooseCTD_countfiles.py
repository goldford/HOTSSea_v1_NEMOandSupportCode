import os

# Set the path to the directory containing the files
path = '/project/6006412/goldford/data/evaluation/nanoose_stn/prepped_pyapnames/'

# Initialize a dictionary to store the count of files by year
file_count_by_year = {}

# Iterate through the files in the directory
for filename in os.listdir(path):
    
    # Split the filename into parts based on the underscore character
    parts = filename.split('_')

    # If the file follows the expected format and has a year as the second part
    if len(parts) == 5 and len(parts[-2]) == 10:
        
        year = int(parts[-2][0:4])
        
        # If the year is not already in the dictionary, add it with a count of 1
        if year not in file_count_by_year:
            file_count_by_year[year] = 1
        # Otherwise, increment the count for that year
        else:
            file_count_by_year[year] += 1

sorted_dict = {k: v for k, v in sorted(file_count_by_year.items())}

# Print the count of files by year
print(sorted_dict)
