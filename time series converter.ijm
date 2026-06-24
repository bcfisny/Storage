//directory containing the images. Use \\ on windows and end the path with \\
dir = "D:\\test macro\\"

// Get a list of all  files in the directory
list = getFileList(dir);

//output directory
dir_out = "D:\\test macro\\output\\"

searchword = "Cy5";

image_num = 0;

// Loop through all files in the directory
for (i = 0; i < list.length; i++) {
	
	
    filePath = "" + dir + "" + list[i];
    
    //skips any folders in the path
    if (!File.isDirectory(filePath)) {
    
    
	    //only opens images that contain a specific word
	    if (indexOf(filePath, searchword) != -1) {
	    	//opens all images
		    open(filePath);
		    
		    image_num += 1;
	    	}
    }
    
}


//combines the images into a stack, saves it and then closes 
run("Images to Stack", "use");

//(channel, Z, Time)
Stack.setDimensions(1, 1, image_num);
run("Bio-Formats Exporter", "save=[" + dir_out + searchword + " full_stack.tif] windowless=true");
run("Close All");



