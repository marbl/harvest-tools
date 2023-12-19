import os,sys,string,subprocess
##harvest utility file for writing VCF file
#ROADMAP
#1. add indel support

helpme = ["-h","-H","-help","-HELP"]

harvest_found = True

#usage/help
if len(sys.argv) == 1 or sys.argv[1] in helpme:
    print("usage: harvest-variants.py <parsnp.hvt> [outfile]")
    sys.exit(1)

#check in path for harvest
try: 
    retcode = subprocess.check_call('harvest',shell=True,stderr=subprocess.PIPE)
except subprocess.CalledProcessError as e: 
    harvest_found = False
except OSError as e:
    harvest_found = False

#if not found, check in HARVESTPATH, otherwise assume currdir
try:
    hpath = os.environ["HARVESTPATH"]
except KeyError:
    hpath = "./"
    if not os.path.exists("./harvest"):
        sys.stderr.write("ERROR: harvest binary couldn't be found, specify location with $HARVESTPATH or place link in current directory")
 
#default output dir, handle stdout/-
outf = "out.vcf"
if len(sys.argv) > 2 and len(sys.argv[2]) >= 1:
    outf = sys.argv[2]

#call harvest
if os.path.exists(sys.argv[1]):
    if harvest_found:
        
        os.system("harvest -q -i %s -V %s"%(sys.argv[1],outf))
    else:
        os.system("%sharvest -q -i %s -V %s"%(hpath,sys.argv[1],outf))
