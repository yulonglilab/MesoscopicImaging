//Fei Deng, 20210525, Add labels
// running time：20220526
// add labels ################################################
// 20210621,增加旋转截图功能
// input parameters ################################################
FrameT = 0.2; // period of each frame/s
BaseFrameNum = 50;  // baseline的帧数
//dirinput = getDirectory("Choose a Directory ");
//diroutput = dirinput;
//print(dirinput);
//print(diroutput);
title = getTitle();
diroutput = getDir("image");
print(diroutput);
getDimensions(width, height, channels, slices, frames);
//print(width +", "+ height +", "+ channels +", "+ slices +", "+ frames);
//print(frames);
// ####### for calcium indicator
// PeakFrameStart = 17; // start of peak response
// PeakFrameEnd = 19; // end of peak response
// ####### for 5-HT3.0
print(title);
TrialNote = substring(title,indexOf(title, "Hz ")+3,indexOf(title, "s_TrialAveRsp.tif")); // 根据文件名提取刺激信息
// print(TrialNote);
PeakFrame = BaseFrameNum+TrialNote/FrameT;
print(PeakFrame);
PeakFrameStart = PeakFrame-1; // start of peak response
print(PeakFrameStart);
PeakFrameEnd = PeakFrameStart+3; // end of peak response
print(PeakFrameEnd);
// frame number of each stage
startS1 = 1; // start of stage1
endS1 = BaseFrameNum; // end of stage1
startS2 = BaseFrameNum+1; // start of stage2
endS2 = PeakFrame; // end of stage2
startS3 = PeakFrame+1; // start of stage3
//endS3 = slices; // end of stage3
endS3 = frames; // end of stage3
// startS4 = 21; // start of stage4
// endS4 = 29; // end of stage4
// startS5 = 30; // start of stage4
// endS5 = frames; // end of stage4

// time of each stage
startS1T = (startS1-BaseFrameNum)*FrameT; // start of stage1
endS1T = (endS1-BaseFrameNum)*FrameT; // end of stage1
startS2T = (startS2-BaseFrameNum)*FrameT; // start of stage2
endS2T = (endS2-BaseFrameNum)*FrameT; // end of stage2
startS3T = (startS3-BaseFrameNum)*FrameT; // start of stage3
endS3T = (endS3-BaseFrameNum)*FrameT; // end of stage3
// startS4T = startS4*FrameT; // start of stage4
// endS4T = endS4*FrameT; // end of stage4
// startS5T = startS5*FrameT; // start of stage4
// endS5T = endS5*FrameT; // end of stage4

// ####### for calcium indicator
// pseudoscaleMin = -0.03; // 伪彩图显示最小值
// pseudoscaleMax = 0.05; // 伪彩图显示最大值
// ####### for 5-HT3.0
// pseudoscaleMin = -0.01; // 伪彩图显示最小值
// pseudoscaleMax = 0.02; // 伪彩图显示最大值
pseudoscaleMin = 0; // 伪彩图显示最小值
// pseudoscaleMax = 2.5; // 伪彩图显示最大值
pseudoscaleMax = 1.6; // 伪彩图显示最大值
// running ################################################
run("Jet");
setMinAndMax(pseudoscaleMin, pseudoscaleMax);

// run("Rotate 90 Degrees Right");
makeRectangle(39, 13, 282, 293);
// makePolygon(38.5,13.3,323.4,13.3,323.4,319.2,38.5,319.2);
run("Crop");
//###### save averaged peak response image  ######
run("Z Project...", "start=" + PeakFrameStart + " stop=" + PeakFrameEnd + " projection=[Average Intensity]");
run("Jet");
setMinAndMax(pseudoscaleMin, pseudoscaleMax);
fullname = diroutput + "/" + substring(title,0,indexOf(title, ".tif")) + "_Frame" + PeakFrameStart + "-" + PeakFrameEnd + "_Peak Ave" + "_Min" + pseudoscaleMin + "Max" +  pseudoscaleMax;
saveAs("Tiff", fullname);
saveAs("png", fullname);

//###### save averaged basal response image  ######
selectWindow(title);
run("Z Project...", "start=" + 1 + " stop=" + BaseFrameNum + " projection=[Average Intensity]");
run("Jet");
setMinAndMax(pseudoscaleMin, pseudoscaleMax);
fullname = diroutput + "/" + substring(title,0,indexOf(title, ".tif")) + "_Frame" + 1 + "-" + BaseFrameNum + "_Basal Ave";
saveAs("Tiff", fullname);
saveAs("png", fullname);

selectWindow(title);
run("RGB Color");
//###### Stage1 ######
setForegroundColor(255, 255, 255);
run("Label...", "format=Text starting=0 interval=" + FrameT + " x=10 y=60 font=25 text=[] range=" + startS1 + "-" + endS1);
run("Label...", "format=00:00 starting=" + startS1T + " interval=" + FrameT + " x=10 y=30 font=25 text=[] range=" + startS1 + "-" + endS1);
//run("Label...", "format=00 starting=" + startS1T + " interval=" + FrameT + " x=10 y=30 font=25 text=[] range=" + startS1 + "-" + endS1);

//###### Stage2 ######
setForegroundColor(255, 0, 0);
run("Label...", "format=Text starting=0 interval=" + FrameT + " x=10 y=60 font=25 text=[Laser] range=" + startS2 + "-" + endS2);
run("Label...", "format=00:00 starting=" + startS2T + " interval=" + FrameT + " x=10 y=30 font=25 text=[] range=" + startS2 + "-" + endS2);
//###### Stage3 ######
setForegroundColor(255, 255, 255);
run("Label...", "format=Text starting=0 interval=" + FrameT + "  x=10 y=60 font=25 text=[] range=" + startS3 + "-" + endS3);
run("Label...", "format=00:00 starting=" + startS3T + " interval=" + FrameT + " x=10 y=30 font=25 text=[] range=" + startS3 + "-" + endS3);

//###### Stage4 ######
// setForegroundColor(0, 255, 0);
// run("Label...", "format=Text starting=0 interval=1 x=10 y=60 font=20 text=[Time window] range=" + startS4 + "-" + endS4);
// run("Label...", "format=00:00 starting=" + startS4T + " interval=" + FrameT + " x=10 y=30 font=20 text=[] range=" + startS4 + "-" + endS4);
//###### Stage5 ######
// setForegroundColor(0, 255, 0);
// run("Label...", "format=Text starting=0 interval=1 x=10 y=60 font=20 text=[Interval] range=" + startS5 + "-" + endS5);
// run("Label...", "format=00:00 starting=" + startS5T + " interval=" + FrameT + " x=10 y=30 font=20 text=[] range=" + startS5 + "-" + endS5);

// save movies ################################################
startSlice = startS1;
//endSlice = slices;
endSlice = frames;
// print(startSlice);
FrameInterNum = 4;
run("Make Substack...", "  slices="+startSlice + "-" + endSlice + "-" + FrameInterNum);
fullname = diroutput + "/" + substring(title,0,indexOf(title, ".tif")) + "_Min" + pseudoscaleMin + "Max" +  pseudoscaleMax + ".avi";
run("AVI... ", "compression=JPEG frame=10 save=fullname");
print("Finished!");
//close("*");