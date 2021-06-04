from cv2 import cv2
import numpy as np
import struct
import os
from bitstring import Bits
for i in range(35):
    file_path = "D:\HEVC\INTRA\ConsoleApplication1\\" + str(i) + ".txt"
    picture = "Mode" + str(i) +".png"
    print(file_path)
    Script     = open(file_path, "r")
    x = 8
    y = 8
    python_img  = np.zeros((x , y))
    for i in range(x):
        for j in range(y):
            python_line       = Script.readline()
            python_img[i][j]  = python_line
#cv2.imshow('Python', python_img)
    cv2.imwrite(picture, python_img)
#cv2.waitKey(0)
#cv2.destroyAllWindows()


