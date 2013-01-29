import distutils.command.sdist
import distutils.command.install

import glob
import os
import tarfile
import shutil

# prevent python from trying to hard link on AFS
delattr(os, "link")

auto_dirs = []

class auto_sdist(distutils.command.sdist.sdist):
	"""
	Custom script for distribution of source
	"""
	global auto_dirs

	def run(self):
		"""
		Main run method.  Steps are as follows:
		1) Configure the directories given
		2) Run make distcheck and ensure it was sucessful
		3) Add the files in the source distribution to the manifest
		4) Run the standard distutils sdist
		"""
		# copy the MANIFEST.in if it exists
		if os.path.exists("MANIFEST.in"):
			shutil.copy2("MANIFEST.in", ".MANIFEST.in_bak")
		
		manifest = file("MANIFEST.in", "a")
		# Print a newline to the manifest, just to make sure we're on a new line
		print >> manifest, " "
		# Of course, add this to the manifest!
		print >> manifest, "include", "autodist.py"
				
		currwd = os.getcwd()
		for d in auto_dirs:
			os.chdir(d)
			os.system("autoreconf -i")
			os.system("./configure")
			os.system("make distcheck")
			tgz_list = glob.glob("*.tar.gz")
			if len(tgz_list) != 1:
				raise Exception("Cannot find a single automake source dist, exiting")
			
			tf = tarfile.open(tgz_list[0],'r')
			for fn in tf.getnames():
				fn_strip = fn.split(os.sep,1)[-1]
				if fn_strip and fn_strip != fn:
					print fn_strip
					print >> manifest, "include", os.path.join(d, fn_strip)
			
			os.chdir(currwd)
		
		manifest.close()
			
		distutils.command.sdist.sdist.run(self)
		
		# Undo the changes we made to the manifest
		if os.path.exists(".MANIFEST.in_bak"):
			shutil.copy2(".MANIFEST.in_bak", "MANIFEST.in")
			os.remove(".MANIFEST.in_bak")
		else:
			os.remove("MANIFEST.in")

class auto_install(distutils.command.install.install):
	global auto_dirs
	
	distutils.command.install.install.user_options.append(("configure=", None, "Additional arguments to pass to the configure script (--prefix and --exec-prefix will be passed automatically)"))
	
	def initialize_options(self):
		distutils.command.install.install.initialize_options(self)
		self.configure = ''
		self.ldprofile = False
	
	def run(self):
		"""
		Main running method.  Algorithm is as follows:
		1) Run standard install
		2) Run configure w/ appropriate options
		3) Run make
		4) Run make install
		"""
		
		distutils.command.install.install.run(self)
		
		currwd = os.getcwd()
		for d in auto_dirs:
			os.chdir(d)
			os.system("./configure --prefix='%s' --exec-prefix='%s' " % (self.prefix, self.exec_prefix) + self.configure)
			os.system("make clean")
			os.system("make")
			os.system("make install")
			os.chdir(currwd)
