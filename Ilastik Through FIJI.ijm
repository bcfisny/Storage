//ensures all operations are performed in the proper order
setBatchMode(true);

//specifies whcih h5 file to use
h5_suffix = "/Data"

//path to ilastik project
ilastik_project = "C:\\Users\\Downloads\\New folder\\test.ilp"

//directory containing the images. Use \\ on windows and end the path with \\
dir = "C:\\Users\\Downloads\\New folder\\h5 folder\\" 

dir_out = "C:\\Users\\Downloads\\New folder\\out put\\"  

// Get a list of all  files in the directory
list = getFileList(dir);

// Loop through all h5 files in the directory
for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], ".h5") || endsWith(list[i], ".hdf5")) {
	
		//save the file name and path
	    filePath =  dir + list[i];
	   
		run("Import HDF5", "select=[" + filePath + "] datasetname=["+ h5_suffix +"] axisorder=tzyxc");
	
		//runs the ilastik pixel classifier
		run("Run Pixel Classification Prediction", "projectfilename=["+ ilastik_project +"] inputimage=["+ list[i] +"] pixelclassificationtype=Probabilities");
	
		//saves the image as a tiff
		saveAs("Tiff", dir_out + list[i]);
	
	    // Close the hdf5 and tiff image after processing
	    run("Close All");
	    
	}
}

//ends batch mode when the code is complete
setBatchMode(false);

