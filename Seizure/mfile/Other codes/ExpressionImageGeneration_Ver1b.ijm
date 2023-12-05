// ExpressionImageGeneration_Ver1
// running time：20220411
//########################################################################################################
//Fei Deng,20210124, mask generation
// roiManager("Select", 1);
title = getTitle();
print(title);
diroutput = getDir("image");

makePolygon(186.1174,28.7,196.357,29.7969,210.6202,33.0883,221.5913,38.5742,232.5631,45.5231,243.5342,55.0312,251.58,64.54,260.3573,75.8772,268.037,87.9459,274.6198,99.6485,279.0088,112.8141,284.494,126.3458,290.346,140.6083,295.4658,152.677,299.4887,165.8426,303.877,183.0311,307.1684,195.4659,310.8259,215.2143,314.8488,233.4997,316.6772,250.6889,315.5803,266.7798,311.5,277.634,304.5,284.634,294.7,286.734,284.2,286.734,268.8,287.434,256.2,287.434,240.8,287.434,224,288.134,205.1,287.434,189.4088,287.2597,172.9511,286.8943,160.3,286.734,143.6946,286.1628,129.7968,285.7974,112.9744,285.7974,100.1742,284.3344,88.4716,283.9683,72.8,281.834,58.8,277.634,51.8,269.234,49,254.534,46.9,241.934,46.9,232.834,47.8772,215.9458,49.3402,207.1685,52.6316,193.2714,55.923,182.3003,57.386,174.2545,59.9459,162.9173,63.9688,150.1171,67.2602,142.8028,70.1855,135.8546,74.2084,124.1513,78.5974,115.0086,82.2542,105.8659,85.1802,98.5516,89.2031,91.6027,93.226,83.5569,101.2711,73.3173,110.4145,61.614,121.3856,50.6429,134.4259,40.334,149.5459,33.4544,160.517,30.5284,172.9511,29.0654);

roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "RoiSet_resized");
fullname = diroutput + "/RoiSet_resized.roi" ;
roiManager("Save", fullname);

run("Make Inverse");
//run("Set...", "value=0");
run("Set...", "value=NaN");
roiManager("Set Color", "white");
roiManager("Set Line Width", 1.5);

makeRectangle(40, 9, 282, 293);
run("Crop");
run("Scale Bar...", "width=1000 height=5 font=14 color=White background=None location=[Upper Right] bold hide overlay");


fullname = diroutput + "/" + substring(title,0,indexOf(title, ".tif")) + "_show" ;
saveAs("Tiff", fullname);
run("Flatten");
run("Copy to System");

close("*");
close("ROI Manager");
print("Finished!");