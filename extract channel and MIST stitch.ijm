//ensures all operations are performed in the proper order
setBatchMode(true);

//use this to change between images where only one channel is needed vs. images where all channels are needed
multichannel = "no"
//multichannel = "yes"

if(multichannel != "yes" && multichannel != "no"){
	print("ERROR. multichannel mode not set properly.");
	exit();
}

//directory containing the raw tiles. Use \\ on windows and end the path with \\
dir = "C:\\Users\\Downloads\\New folder\\";

// Get a list of all  files in the directory
list = getFileList(dir);



//
//inputs for multi vs. singlechannel mode
//

if(multichannel == "no"){

	//directory to save new single channel images
	dir_out = "C:\\Users\\Downloads\\New folder\\Single_Channel\\";
	
	// Check if the strings are identical
	if (dir == dir_out) {
	    exit("Error: The input and output path are the same. This was cause unnecessary files to be duplicated.");
	}
	
	//Adds a prefix to the images with only one channel
	prefix = "BLUE";
	
	//choose which RGB (1,2 or 3) channel to extract
	channel = 3;

	// Loop through all files in the directory
	for (i = 0; i < list.length; i++) {
		
		// Build full path
	    filePath =  dir + list[i];
	    
	    //skips any folders in the path
    	if (!File.isDirectory(filePath)) {
    		
		    run("Bio-Formats Importer", "open=[" + filePath + "]");
		
		    // Get the image title (needed for channel extraction)
		    title = getTitle();
		
			// Split channels
			run("Split Channels");
		
			// Channel window names created by Split Channels
			if (channel == 1)
		    		selectWindow("C1-" + title);
			else if (channel == 2)
		    		selectWindow("C2-" + title);
			else if (channel == 3)
		    		selectWindow("C3-" + title);
			else
		    		exit("Invalid channel number");
		    
		    
		    saveAs("Tiff", dir_out + prefix + list[i]);
		    close("*");
    	}
	    
	}

}


//uses the same directory as the raw tiles if all channels are needed
else if(multichannel == "yes"){
	dir_out = dir;
	prefix = "";
}

//
//inputs for stitching
//

//directory to save stitched image output
dir_final = "E:\\data\\Ben\\RGB_Test\\Output2\\";

//use p for numbers and add .tif to the end (try .tif.tif if the filename ends in .tif)
filename = prefix + "AD3C-2-20X_{ppppp}_CH4.tif";

//number of tiles in each direction
horizontal_width = 18;
vertical_width = 14;

//pixel dimensions for Keyence
/*
 * objective	width/height (um)
 * 2x			3.77
 * 4x			1.88
 * 10x			0.755
 * 20x			0.377
 * 40x			0.18872
 * 100x			0.07549
*/

width_px = 0.377;
height_px = width_px;


setBatchMode(false);

//Stitches the images together and puts the final image in one folder
run("MIST",
    "gridwidth=" + horizontal_width + " gridheight="+ vertical_width +" starttile=1 " +
    "imagedir=[" + dir_out + "] " +
    "filenamepattern=[" + filename + "] " +
    "filenamepatterntype=SEQUENTIAL " +
    "gridorigin=UL assemblefrommetadata=false assemblenooverlap=false " +
    "globalpositionsfile=[] numberingpattern=HORIZONTALCONTINUOUS " +
    "startrow=0 startcol=0 " +
    "extentwidth=" + horizontal_width + " " +
    "extentheight=" + vertical_width + " " +
    "timeslices=0 istimeslicesenabled=false " +
    "outputpath=[" + dir_final + "] " +
    "displaystitching=false outputfullimage=true outputmeta=true " +
    "outputimgpyramid=true blendingmode=OVERLAY blendingalpha=NaN " +
    "compressionmode=ZLIB outfileprefix=stitched- " +
    "unit=MICROMETER unitx=" + width_px + " unity=" + height_px + " " +
    "programtype=AUTO numcputhreads=32 loadfftwplan=true " +
    "savefftwplan=true fftwplantype=MEASURE " +
    "fftwlibraryname=libfftw3 fftwlibraryfilename=libfftw3.dll " +
    "planpath=[C:\\Users\\Downloads\\Fiji.app\\lib\\fftw\\fftPlans] " +
    "fftwlibrarypath=[C:\\Users\\Downloads\\Fiji.app\\lib\\fftw] " +
    "stagerepeatability=0 horizontaloverlap=NaN verticaloverlap=NaN " +
    "numfftpeaks=0 overlapuncertainty=NaN " +
    "isusedoubleprecision=false isusebioformats=false " +
    "issuppressmodelwarningdialog=false isenablecudaexceptions=false " +
    "translationrefinementmethod=SINGLE_HILL_CLIMB " +
    "numtranslationrefinementstartpoints=16 " +
    "headless=false loglevel=MANDATORY debuglevel=NONE"
);