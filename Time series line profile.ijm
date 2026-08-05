//ensures all operations are performed in the proper order
setBatchMode(true);

id = getImageID();
getDimensions(width, height, channels, slices, frames);

//Use \\ on windows and end the path with \\
dir_out = "C:\\Users\\Downloads\\New folder\\out put\\";




//loop through each frame
for (frame = 0; frame < frames; frame++) {

    selectImage(id);
    Stack.setPosition(1, 1, frame+1);

	//gets the intensity profile of the line
    profile = getProfile();

	//adds the T slice number to the data
    setResult("Frame", frame, frame+1);
   

	//loops through every value in the intensity profile so they can all be saved per frame
    for (p = 0; p < profile.length; p++) {
        setResult("Intensity "+p, frame, profile[p]);
    }
}

updateResults();
saveAs("Results", dir_out + "Intensity Results.csv");




//get the length of the line
run("Measure");
lineLength = getResult("Length", nResults-1);

//get the number of points on the line
profile = getProfile();
n = profile.length;

//line step does not include the first point 
line_step = lineLength/(n-1)

run("Clear Results");

//used to properly align the columns of the distance measurements to the intensity values
setResult("", 0 , "")

//loops through all intensity values and saves the distance measurements
for (i = 0; i < profile.length; i++) {
    setResult("Distance step " + i, 0, line_step*i);
}

updateResults();
saveAs("Results", dir_out + "Distance Results.csv");




dist = File.openAsString(dir_out + "Distance Results.csv");
intensity = File.openAsString(dir_out + "Intensity Results.csv");

title = getTitle();

//combines the distance and intensity data together
File.saveString(dist + intensity, dir_out + title + " Combined Results.csv");

//closes any open windows
run("Close");

//ends batch mode when the code is complete
setBatchMode(false);