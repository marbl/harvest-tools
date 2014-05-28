#!python
import os,sys,string,subprocess

##harvest utility file for multi-fasta file of SNP loci

##ROADMAP
#1. add option to output indels
#2. add option to filter columns?

helpme = ["-h","-H","-help","-HELP"]
harvest_found = True

#print help/usage
if len(sys.argv) == 1 or sys.argv[1] in helpme:
    print "usage: harvest-snps.py <parsnp.hvt> [outfile]"
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
outf = "out.vcf"

if len(sys.argv) > 2 and len(sys.argv[2]) >= 1:
    outf = sys.argv[2]

if os.path.exists(sys.argv[1]):# or sys.argv[1] == "stdout":
    if harvest_found:
        
        os.system("harvest -q -i %s -S %s"%(sys.argv[1],outf))
    else:
        os.system("%sharvest -q -i %s -S %s"%(hpath,sys.argv[1],outf))
