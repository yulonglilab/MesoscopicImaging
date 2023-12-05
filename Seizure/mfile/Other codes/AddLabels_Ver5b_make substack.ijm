//AddLabels_Ver5b_make substack, Fei Deng
// running time：20220501
// add labels ################################################
// 20220206,使用自定义函数简化代码
// 20220501,导出视频前先转换为RGB格式
// input parameters ################################################
FrameT = 1; // period of each frame/s
title = getTitle();
getDimensions(width, height, channels, slices, frames);
//print(width +", "+ height +", "+ channels +", "+ slices +", "+ frames);
dirinput = getDir("image");
diroutput = dirinput + "/Waves"; // 截取视频存储文件夹
print(diroutput);
File.makeDirectory(diroutput);
//print(frames);
// ####### for calcium indicator
// PeakFrameStart = 17; // start of peak response
// PeakFrameEnd = 19; // end of peak response
// ####### for 5-HT3.0
print(title);
FrameInterNum = 1;
// ####### for jRGECO1a sensor
//pseudoscaleMin = -0.1; // 伪彩图显示最小值
//pseudoscaleMax = 0.9; // 伪彩图显示最大值
// ####### for 5-HT3.0mut sensor
//pseudoscaleMin = -0.2; // 伪彩图显示最小值
//pseudoscaleMax = 0.5; // 伪彩图显示最大值

// ####### for 5-HT3.0 sensor
pseudoscaleMin = -0.3; // 伪彩图显示最小值
pseudoscaleMax = 2; // 伪彩图显示最大值
// ####### for r5-HT2.0 sensor
//pseudoscaleMin = -0.2; // 伪彩图显示最小值
//pseudoscaleMax = 0.4; // 伪彩图显示最大值
// ####### for eCB2.0 sensor
//pseudoscaleMin = -0.2; // 伪彩图显示最小值
//pseudoscaleMax = 4; // 伪彩图显示最大值
run("Jet");
setMinAndMax(pseudoscaleMin, pseudoscaleMax);

//run("Rotate 90 Degrees Right");
makeRectangle(37, 16, 282, 293);
run("Crop");
//###### save averaged basal response image  ######
//  selectWindow(title);
// run("Z Project...", "start=" + 1 + " stop=" + BaseFrameNum + " projection=[Average Intensity]");
// run("Jet");
// setMinAndMax(pseudoscaleMin, pseudoscaleMax);
// fullname = diroutput + "/" + substring(title,0,indexOf(title, ".tif")) + "_Frame" + 1 + "-" + BaseFrameNum + "_Basal Ave";
// saveAs("Tiff", fullname);
// saveAs("png", fullname);

// 手动输入每个block的起始帧数
BlockListS = newArray(19401,0);  // 由于newArray只有一个元素时，BlockListS.length和Array的索引结果有问题，故在末尾加一个0
BlockListE = newArray(19600,0);
// BlockListS = newArray(7461,7861,7961,8261,8511,8790,13411,13881,14111,14241,14361,14511,14631);
// BlockListE = newArray(7670,7940,8046,8360,8660,8959,13630,13990,14170,14380,14460,14600,16245);

if (BlockListS[BlockListS.length-1] == 0) {
	BlockNum = 1; } 	
else{
	BlockNum = BlockListS.length;
}
print("Block number: "+BlockNum);

//###### Generate blocks  ################################################################################################################
for (i=0; i<BlockNum; i++) {
  blockNo = i+1;
  print("Generating block "+blockNo);
  startSlice = BlockListS[i]; // start of stage1
  endSlice = BlockListE[i]; // end of stage1
  print("Start fame from: "+startSlice+", End fame from: "+endSlice);
  MakeSubstack(title, startSlice, endSlice, FrameInterNum,pseudoscaleMin,pseudoscaleMax);
}
//######################################################################################################################################################

print("Finished!");
//close("*");

//######################################################################################################################################################
function MakeSubstack(title, startSlice, endSlice, FrameInterNum, pseudoscaleMin, pseudoscaleMax) {
  selectWindow(title);
  run("Make Substack...", "  slices="+startSlice + "-" + endSlice + "-" + FrameInterNum);
  fullname = diroutput + "/" + substring(title,0,indexOf(title, ".tif")) + "_Frame" + startSlice + "-" + endSlice +"_Min" + pseudoscaleMin + "Max" +  pseudoscaleMax + ".tif";
  saveAs("Tiff", fullname);
  titleSub = getTitle();
  //###### Stage n ######
  run("RGB Color");
  setForegroundColor(255, 255, 255);
  startSliceSub = startSlice-startSlice+1;
  endSliceSub = endSlice-startSlice+1;
  //run("Label...", "format=Text starting=0 interval=" + FrameT + " x=10 y=60 font=25 text=[] range=" + startS1 + "-" + endS1);
  run("Label...", "format=0 starting=" + 0 + " interval=" + FrameT + " x=10 y=30 font=25 text=[s] range=" + startSliceSub + "-" + endSliceSub);

  // save movies ################################################
  // FrameInterNum = 1;
  // run("Make Substack...", "  slices="+startSlice + "-" + endSlice + "-" + FrameInterNum);
  fullname = diroutput + "/" + substring(title,0,indexOf(title, ".tif")) + "_Frame" + startSlice + "-" + endSlice +"_Min" + pseudoscaleMin + "Max" +  pseudoscaleMax + ".avi";
  run("AVI... ", "compression=JPEG frame=10 save=fullname");
}