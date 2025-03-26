#!/usr/bin/python3

#Process output file from Slattice into something readable

import sys

def main(filename):
  # define output file name
  outputfile = filename + ".simple"
  data = 0

  # open the output file
  with open(outputfile,"w+") as out_short:
    # open the original slattice output
    with open(filename) as out_full:
      # going line by line
      for line in out_full:
        # if the line starts with "[1]", it will contain results
        if line[0:3] == "[1]":
          # split the line
          elements = line.split()
          # ckeck if the third element is a selection coefficient and write it out
          if elements[3] == "s:":
            data += 1
            newline = elements[4][:-1] + "\n"
            out_short.write(newline)
          # check if the third element is a confidence interval ant write it out
          elif elements[3] == "CI:":
            data += 2
            newline = elements[4] + "\n"
            out_short.write(newline)
            newline = elements[5][:-1] + "\n"
            out_short.write(newline)
      # if there was no result after going through all the lines, output NAs
      if data == 0:
        newline = "NA\nNA\nNA\n"
        out_short.write(newline)
      # if there was only a selection coefficient, but no CI, fill with 2 trailing NAs
      elif data == 1:
        newline = "NA\nNA\n"
        out_short.write(newline)


if __name__ == "__main__":
  arg1 = sys.argv[1]
  main(arg1)
