
import os
import pickle
in_path = "/Users/weeder/PycharmProjects/pepsickle-paper/data/training_sets"

files = os.listdir(in_path)


for file in files:
    if file.endswith(".pickle"):
        in_file = in_path + "/" + file
        dat = pickle.load(open(in_file, "rb"))

        print(file)
        if file.startswith("all_mammal_20S"):
            pos = len(dat['proteasome']['positives'])
            neg = len(dat['proteasome']['negatives'])

            print("Positives in set: ", pos)
            print("Negatives in set: ", neg)
        elif file.startswith("all_mammal_epitope"):
            pos = len(dat['epitope']['positives'])
            neg = len(dat['epitope']['negatives'])

            print("Positives in set: ", pos)
            print("Negatives in set: ", neg)

        else:
            print("doesn't fit regex")

