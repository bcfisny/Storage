/* Script designed to detect cells in H&E images
 * The user creates a manual annotation over the area of interest and classifies it as "Circle"
 * Then cells are detected over the entire image and classified based on their eosin intensity
 * Cells located in the background of the image, outside of the eosin, are also included inside the "Positive" annotation
 *The "Negative" annotation is split to try and find the area of individual eosin branches 
 */


setImageType('BRIGHTFIELD_H_E')
setColorDeconvolutionStains('{"Name" : "H&E default", "Stain 1" : "Hematoxylin", "Values 1" : "0.6511078257574492 0.7011930431234068 0.29049426072255424", "Stain 2" : "Eosin", "Values 2" : "0.21589893562087106 0.8011960501132093 0.5580972485873467", "Background" : " 255 255 255"}');
setPixelSizeMicrons(0.377, 0.377)

//removes any annotations that are not labeled. This avoids duplicates
selectObjects { it.isAnnotation() && it.getPathClass() == null}
remove_null = getSelectedObjects()
removeObjects(remove_null, true)
resetSelection();

//removes any annotations that are not labeled. This avoids duplicates
selectObjectsByClassification("Positive", "Negative")
remove = getSelectedObjects()
removeObjects(remove, true)
resetSelection();


selectObjectsByClassification("Circle")
makeInverseAnnotation()

createAnnotationsFromPixelClassifier("Background threshold", 150.0, 150.0)

selectObjectsByClassification("Circle","Positive")
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImageBrightfield":"Hematoxylin OD","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":9.8,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":0.5,"minAreaMicrons":10.0,"maxAreaMicrons":75.0,"threshold":0.1,"maxBackground":2.0,"watershedPostProcess":true,"cellExpansionMicrons":2.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}')
runObjectClassifier("eosin classifier")

selectObjectsByClassification("Negative")
runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin', '{}')


//extracts the name of the image so it can be used when saving the output statistics 
//def imageName = GeneralTools.getNameWithoutExtension(getCurrentImageData().getServer().getMetadata().getName())

//saves the statistics to a csv (must use forward slashes for the path on windows)
//saveDetectionMeasurements('/path'+imageName+'.csv')