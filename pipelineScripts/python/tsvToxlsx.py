#!/bin/python
import argparse
from openpyxl import Workbook
import csv
import os

parser = argparse.ArgumentParser(description="transform a tab separated text to an Excel-ready file", epilog="Ouputfile will be: textfile.xlsx")
parser.add_argument('-i', '--inputfile', help='[REQUIRED] A tab separated text', dest='inputfile', action='store', required=True)
args = parser.parse_args()
inputfile = args.inputfile

prefix = inputfile.split('.')[0]
outputfile = prefix + '.xlsx'

os.system('sed \"s/,/;/g\" ' + inputfile + '> _temp')
os.system('sed -i \"s/\t/,/g\" _temp')


wb = Workbook()
ws = wb.active
with open('_temp', 'r') as f:
    for row in csv.reader(f):
        ws.append(row)
wb.save(outputfile)


os.system('rm -f _temp')
