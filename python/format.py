#!/usr/bin/env python
# Use python to convert the dictionary format data into rows of data.

from pandas import Series,DataFrame

myfile = open("/mnt/harddrive/light/a.txt")
myfile.read()
str1 = myfile.read()
myfile.seek(0)
str1 = myfile.read()
str1.split()
data = str1
frame1 = DataFrame(eval(data))

frame1.to_csv("/mnt/harddrive/light/b.txt", '\t', index = False)

