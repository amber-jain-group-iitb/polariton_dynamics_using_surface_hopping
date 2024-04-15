import os
import time
import shutil

def copying_files(i):
	os.mkdir('output'+str(i))
	shutil.copy("AFSSH.f90",'output'+str(i))
	shutil.copy("AFSSH.inp",'output'+str(i))
	shutil.copy("AFSSH.o",'output'+str(i))
	#shutil.copy("20_knot_points_w.txt",'output'+str(i))
	#shutil.copy("20_knot_points_x.txt",'output'+str(i))
	#shutil.copy("knot_points_x.inp",'output'+str(i))
	#shutil.copy("knot_points_w.inp",'output'+str(i))
	shutil.copy("mod_afssh.f90",'output'+str(i))
	shutil.copy("mod_afssh.mod",'output'+str(i))
	shutil.copy("mod_afssh.o",'output'+str(i))
	shutil.copy("sub.sh",'output'+str(i))
	shutil.copy("aout",'output'+str(i))
	os.chdir('output'+str(i))
	file1=open("ifolder.inp","w")
	file1.write(str(i))
	file1.close()
	os.system('qsub sub.sh')
	#os.system('qsub -l nodes=node6 sub.sh')
	os.chdir('..')

#os.system("rm -r output*")
#for i in range(1,41):
#	os.system("rm -r output"+str(i))

for i in range(41,81):
	copying_files(i)
#	os.system("rm -r output"+str(i))


#1-40 4000K temp
#41-80 5000K temp
