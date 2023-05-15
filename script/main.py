import argparse
import os

import create_csv
import find_anomaly

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--csv", dest ="csv",default="", help="true/false run csv")
parser.add_argument("-b", "--both", dest ="both", help="true/false run csv and find")

# The parameters for running CSV
parser.add_argument("-exp", "--experiment", dest = "experiment", default = "", help="The path of the patients")
parser.add_argument("-s", "--save",dest ="save", help="The path where the CSV will be saved")

# The parameters for running FIND ANOMALY
parser.add_argument("-df", "--df", dest="df", help="The path of the csv")
parser.add_argument("-num", "--number", dest="number", help="patient number")
parser.add_argument("-pro", "--protein", dest="protein", help="protein name")
parser.add_argument("-mod", "--modelname", dest="modelname", help="XGBoost/DecisionTree")
parser.add_argument("-cell", "--celllabel", dest="celllabel", help="The path of the images cell label")
parser.add_argument("-std",  dest="std", help="")


args = parser.parse_args()
if __name__ == "__main__":
    if args.csv == "True":
         print("run csv")
         print(f'path : {args.save}')

         result_df = create_csv.main(args.experiment, args.save)
         print("created cellTable successfully")

         if args.both == "True":
             print("run find anomaly")
             print(f'df : {result_df}')
             find_anomaly.main(args.both, result_df, args.number, args.protein, args.modelname, args.save)
             print("done!")
         else:
             print("Done!")
    else:
         print("run find anomaly")
         # find_anomaly.main(args.both, args.df, args.number, args.protein, args.modelname, args.celllabel)
         print("done find anomaly")

