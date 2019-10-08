#!/usr/bin/env python3

import os
import sys


def runConvert(fileName):
    if ".raw" in fileName:
        outFile = species[fileName.strip(".raw")] + ".mgf"
        command = "docker run --rm -e WINEDEBUG=-all -v /data1/zhoulab/wukai/data/cotton/data/:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert --mgf /data/other_raw/%s --outfile /data/mgf/%s -o /data/mgf/ ; chown wukai:zhoulab /data/mgf/%s"%(fileName, outFile, outFile)
        os.system(command)

def main():
    if len(sys.argv) != 1:
        sys.stdout.write("Usage:Programe \n")
        sys.exit(1)

    rawFile_path = "data/other_raw/"
    rawFile_list = os.listdir(rawFile_path)

    for r in rawFile_list:
        runConvert(r)  

if __name__ == "__main__":
    main()
