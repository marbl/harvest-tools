#!python3
import os,sys,string,subprocess

##harvest utility file for XMFA file(s) of core genome alignments

##ROADMAP
#1. add option to create multiple files (one per lcb)?
#2. add option to filter on lcb length?
#3. add option to output in MAF format?

helpme = ["-h","-H","-help","-HELP"]
harvest_found = True

#print help/usage
if len(sys.argv) == 1 or sys.argv[1] in helpme:
    print("usage: harvest-alignments.py <parsnp.hvt> [outfile]")
    sys.exit(1)

#check in path for harvest
try: 
    retcode = subprocess.check_call('harvest',shell=True,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e: 
    harvest_found = False
except OSError as e:
    harvest_found = False
try:
    hpath = os.environ["HARVESTPATH"]
except KeyError:
    hpath = "./"
    if not os.path.exists("./harvest"):
        sys.stderr.write("ERROR: harvest binary couldn't be found, specify location with $HARVESTPATH or place link in current directory")

#deafault output file
outf = "out.xmfa"

if len(sys.argv) > 2 and len(sys.argv[2]) >= 1:
    outf = sys.argv[2]

if os.path.exists(sys.argv[1]):# or sys.argv[1] == "stdout":
    if harvest_found:
        
        os.system("harvest -q -i %s -X %s"%(sys.argv[1],outf))
    else:
        os.system("%sharvest -q -i %s -X %s"%(hpath,sys.argv[1],outf))
