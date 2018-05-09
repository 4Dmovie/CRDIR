import rawpy
import imageio
import cv2
import numpy as np
import subprocess

path = 'C:\CRDIR\issueeeeeeee10.02\iss031e152'

startingstring = input('Please type in the starting file name:')
endingstring = input('Please type in the ending file name:')
ImageBit = input('Image bits:')
ImageRows = input('number of rows per image:')
ImageCols = input('number of columns per image:')
text_file = open("H:/RIT intern/codes/test/data.txt", "w")
text_file.write(startingstring + '\n')
text_file.write(endingstring + '\n')
text_file.write(int(ImageBit) + '\n')
text_file.write(int(ImageRows) + '\n')
text_file.write(int(ImageCols) + '\n')
text_file.close()
#np.savetxt('C:\CRDIR\issueeeeeeee10.02\data.txt', )

for x in range(51, 55):
   s = path + str(int(x / 100)) + str(int((x % 100) / 10)) + str(int(x % 10)) + ".nef"
   print (s)
   rawim = rawpy.imread(s)
   temp = rawim.raw_image
   temp = temp.astype(int)
   print(np.shape(temp))
   #np.savetxt(, temp)
   np.savetxt(s + '.txt', temp, fmt='%d')

subprocess.call(["H:/RIT intern/codes/test/Debug/test.exe"])
  # tmp = subprocess.call("./a.out")
   #print
   #"printing result"
   #print
   #tmp
