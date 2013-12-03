import os,sys,string
#harvest utilities INSTALL script
user_home = os.environ["HOME"]
print "<<Welcome to harVest utility script install>>"

#check for python version                                                                                                                                                                         
if (sys.version_info[0] < 2) or (sys.version_info[0] == 2 and sys.version_info[1] < 6):

  print "Python version is %s. metAMOS requires at least 2.6"%(sys.version)
  sys.exit(1)

#complete shebang
os.system("python setup.py install_scripts --install-dir=`pwd` build")

scripts = ["harvest-alignments.py","harvest-phylogeny.py","harvest-backbone.py","harvest-variants.py","harvest-snps.py"]
#copy to currdir
files = os.listdir(".")
for script in scripts:
    os.system("mv %s %s"%(script,script.replace(".py","")))

