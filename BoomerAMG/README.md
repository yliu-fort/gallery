# Hypre2.0 must be compiled as shared library  

### Set environment variables  
	HYPRE_SRC = <Hypre-root-path> >> ~/.bashrc  
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HYPRE_SRC/hypre/lib >> ~/.bashrc

### compile  
	wclean && wmake

### add entry into controlDict to link compiled library  
	libs  
	{  
		"libBoomerAMG.so"	  
	}  

###
### Copyrights 2019 Yuxuan Liu.  
This is an openFOAM extension.  
The goal of this extension is to provide a interface for 3rd-party computing libraries 
like HYPRE.  
This is something like a prototype by now.  
Does not support cyclic boundary currently.  

