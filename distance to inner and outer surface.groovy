//necessary import statements for the interseciton of the two annotations
import qupath.lib.gui.commands.Commands;
import qupath.lib.roi.RoiTools.CombineOp;

//default image values
setImageType('BRIGHTFIELD_H_DAB');
setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049", "Stain 2" : "DAB", "Values 2" : "0.26917 0.56824 0.77759", "Background" : " 255 255 255"}');
resetSelection();

//first selects manual annotations (if they exist) and then fills in red pixels inside those annotations
selectAnnotations();
//Locates red pixels in the image 
createAnnotationsFromPixelClassifier("Red Pixel Classifier", 0.0, 0.0)
selectObjectsByClassification("Not-Red")
oth = getSelectedObjects()
removeObjects(oth, true)
resetSelection();

//Segments the surface
createAnnotationsFromPixelClassifier("Surface Downsampled", 100000.0, 0.0)
//createAnnotationsFromPixelClassifier("K NN D4", 100000.0, 0.0)
//createAnnotationsFromPixelClassifier("Surface", 100000.0, 0.0)
resetSelection();


//Removes small objects and the outer surface from the positive channel so it only lines the inner surface
selectObjectsByClassification("Outside")
runPlugin('qupath.lib.plugins.objects.RefineAnnotationsPlugin', '{"minFragmentSizePixels":100000,"maxHoleSizePixels":1.0E22}')

//selects both the inside and outside channels
out = getAnnotationObjects().findAll{p -> p.getPathClass() == getPathClass("Outside")}[0]
ins = getAnnotationObjects().findAll{p -> p.getPathClass() == getPathClass("Inside")}[0]

//Select the outside channel, then the inside, then take the intersection
//This creates an inside inner semgentation and positive outer segmentation
getCurrentHierarchy().getSelectionModel().setSelectedObject(out, true);
getCurrentHierarchy().getSelectionModel().setSelectedObject(ins, true);
Commands.combineSelectedAnnotations(getCurrentImageData(), CombineOp.INTERSECT);

//Duplicates the positive image since it was removed in the intersection
addObject(out)

//Creates superpixels for the red pixels
selectObjectsByClassification("Red Pixels")
//change the value of "spacingPixels":50.0 to something like "spacingPixels":15.0
runPlugin('qupath.imagej.superpixels.SLICSuperpixelsPlugin', '{"sigmaPixels":5.0,"spacingPixels":50.0,"maxIterations":10,"regularization":0.25,"adaptRegularization":false,"useDeconvolved":false}')

//finds the distance from the superpixels to the outer and inner annotations
detectionToAnnotationDistancesSigned(true)
//Make sure you use forward slashes for the path
saveDetectionMeasurements('C:/Users/bfisler/.openjfx/Documents/QuPath Training')
saveAnnotationMeasurements('C:/Users/bfisler/.openjfx/Documents/QuPath Training')