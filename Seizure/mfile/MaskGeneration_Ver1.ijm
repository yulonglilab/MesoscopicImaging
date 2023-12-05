// MaskGeneration_Ver1
//########################################################################################################
//Fei Deng,20220326, mask generation for images import from min image
// running time：20220503
//#### Input parameters #######################################################################################


//#### Processing #######################################################################################
diroutput = getDir("image");
print(diroutput);
getDimensions(width, height, channels, slices, frames);
print(slices);
// adjust ROI manually
roiManager("Save", diroutput + "/RoiSet.zip");
roiManager("Select", 1);
run("Make Inverse");
//run("Set...", "value=0");
run("Set...", "value=NaN");
run("Trainable Weka Segmentation");
selectWindow("Trainable Weka Segmentation v3.3.2");


run("Set...", "value=0");
run("Set...", "value=1");
       
setThreshold(1, 255);
run("Convert to Mask");
run("Grays");
//run("Save");
saveAs("Tiff", diroutput + "/WholeCortexMask.tif");

close("*");
close("ROI Manager");
print("Finished!");