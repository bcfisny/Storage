//ensures all operations are performed in the proper order
setBatchMode(true);

//directory containing the images. Use \\ on windows
dir = "C:\\Users\\Downloads\\New folder\\";

// Get a list of all  files in the directory
list = getFileList(dir);

//File type to convert to and suffix 
newsuffix = ".h5";

h5_ending = "/Data";

//output directory
dir_out = "C:\\Users\\Downloads\\New folder\\out put\\";


// Loop through all h5 files in the directory
for (i = 0; i < list.length; i++) {
	
	
	// Build full path
    filePath =  dir + list[i];
    
    //skips any folders in the path
    if (!File.isDirectory(filePath)) {
    
	    run("Bio-Formats Importer", "open=[" + filePath + "]");
	    
	    run("Export HDF5", "exportpath=[" + dir_out + list[i] + newsuffix + "] " + "datasetname=[" + h5_ending + "] " + "compressionlevel=0");
	    close("*");
    }
    
}

//ends batch mode when the code is complete
setBatchMode(false);


