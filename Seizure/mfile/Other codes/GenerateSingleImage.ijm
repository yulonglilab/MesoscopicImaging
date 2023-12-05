//GenerateSingleImage, Fei Deng
// running time：20220627
// ################################################
// 20220627,批量生成单帧图像
// input parameters ################################################
FrameT = 1; // period of each frame/s
title = getTitle();
print(title);
getDimensions(width, height, channels, slices, frames);
//print(width +", "+ height +", "+ channels +", "+ slices +", "+ frames);
dirinput = getDir("image");
diroutput = dirinput + "/Examples"; // 截取视频存储文件夹
print(diroutput);
File.makeDirectory(diroutput);
// ####### for 5-HT3.0 sensor
pseudoscaleMin = -0.1; // 伪彩图显示最小值
pseudoscaleMax = 0.4; // 伪彩图显示最大值
// ####### for r5-HT2.0 sensor
//pseudoscaleMin = -0.2; // 伪彩图显示最小值
//pseudoscaleMax = 0.4; // 伪彩图显示最大值
// ####### for eCB2.0 sensor
//pseudoscaleMin = -0.2; // 伪彩图显示最小值
//pseudoscaleMax = 4; // 伪彩图显示最大值
run("Jet");
setMinAndMax(pseudoscaleMin, pseudoscaleMax);

//run("Rotate 90 Degrees Right");
//makeRectangle(44, 18, 282, 293);
//run("Crop");


// BlockListS = newArray(14115,14266,14396,14468,14792,15019,15370,15820,17173,17216,17384);  // 由于newArray只有一个元素时，BlockListS.length和Array的索引结果有问题，故在末尾加一个0
BlockListS = newArray(1783,2004,2148,2277,2359,2375,6550,7960,9555);  // 由于newArray只有一个元素时，BlockListS.length和Array的索引结果有问题，故在末尾加一个0
// REM,Wake,NREM,Wake,NREM,NREM,Wake,NREM,REM
if (BlockListS[BlockListS.length-1] == 0) {
	BlockNum = 1; } 	
else{
	BlockNum = BlockListS.length;
}
print("Crop number: "+BlockNum);

for (i=0; i<BlockNum; i++) {
  blockNo = i+1;
  print("Generating image "+blockNo);
  startSlice = BlockListS[i]; // start of stage1
  print(startSlice);
  //  startSlice = 14115;
  setSlice(startSlice);
  run("Duplicate...", "use");
  print("Save fame : "+startSlice);
  fullname = diroutput + "/" + substring(title,0,indexOf(title, ".tif")) + "_Frame" + startSlice;
  saveAs("Tiff", fullname);
  saveAs("png", fullname);
  close();
}

print("Finished!");
//close("*");
