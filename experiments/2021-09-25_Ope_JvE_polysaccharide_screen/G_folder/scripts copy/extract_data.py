#IF MEASURED WITH OD540 instead of 600, CHANGE LiNE 54

import argparse
import pandas as pd

def parse_args():
    # Build out parser
    parser = argparse.ArgumentParser()
    
    # Add argument for the raw kinetic OD600 data
    parser.add_argument("-i", "--input", 
                        help="txt file with kinetic data from OD600 plate", 
                        type=str, required=True)

    # Add argument for the csv file to write out extracted data to
    parser.add_argument("-o", "--output", 
                        help="csv file to write extracted data from", 
                        type=str, required=True)

    args = parser.parse_args()
    return args


def pull_data(OD600_file):
    # Trigger flag to start storing data
    IS_OD600_DATA = False

    # Trigger flag to stop storing data
    PAST_OD600_DATA = False

    # Stores dicts to convert into dataframe
    df_rows = []

    with open(OD600_file, mode='r', encoding="ISO-8859-1") as f:
        for line in f:
            # until OD600 data is reached, check to see if it is reached
            # then change flag
            if not IS_OD600_DATA:
                if line.strip().startswith("Time	0:00:00"):
                    # Pull out the header
                    header = line.strip().split("\t")[1:]
                    data = pd.DataFrame(columns=header)
                    data.columns.name = 'time'

                    # change flag
                    IS_OD600_DATA = True

            # storable data
            elif line.strip() != "":
                vals = line.strip().split("\t")
                col_name = vals[0]
                
                
                if col_name == "T° 600":
                    col_name = "temp"
                
                vals = vals[1:]
                
                ddict = {}
                for col, val in zip(header, vals):
                    ddict[col] = val
                
                vals = pd.Series(ddict, name = col_name)
                data = data.append(vals)

            # stored all data already
            else:
                break

    data = data.transpose()
    data.index.name = "time"
    return data
        
        
def tidy_data(data):

    # drop NA rows
    data = data.dropna(axis=0, how="all")

    # tidy the data
    data = data.melt(id_vars = ["temp"], var_name = "well", value_name = "OD600", ignore_index = False)
    
    # clean up wells
    data.well = data.well.str[0] + data.well.str[1:].str.zfill(2)
    
    return data

def convert_times(row):
    timestamp = row.name
    
    timestamp = timestamp.split(":")
    
    return(60*int(timestamp[0]) + int(timestamp[1]))


if __name__ == "__main__":
    args = parse_args()

    # Pull out the OD600 file stored as input
    OD600_file = args.input

    # process OD600_file and store the resulting df
    data = pull_data(OD600_file)

    # clean up the data
    data = tidy_data(data)

    # convert the time index to minutes since start
    data.index = data.apply(convert_times, axis=1)
    data.index.name = "time_mins"

    # write it out to out file
    data.to_csv(args.output)



