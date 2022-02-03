/*	210719_BIL_MonicaVandenberg_Neurotransmitters_v2
	***********************
	
	ImageJ/Fiji macro set for to detection of neurotransmitter markers
	Buidling blocks used from CellBlocks.ijm written by Winnok H. De Vos - winnok.devos@uantwerpen.Be
	
	Author: 			Marlies Verschuuren - marlies.verschuuren@uantwerpen.be
	Date Created: 		2021 - 07 - 19
	Date Last Modified:	2022 - 01 - 20 
*/

/*	
	Plugins needed (Help > Update > Manage update sites):
 	
 	Bio-Formats: https://sites.imagej.net/Bio-Formats/
	ImageScience: https://sites.imagej.net/ImageScience/
*/

/*	
 	VERSION SUPPORT: 
	V1:	- Spot detection in different channels
		- Max projection of substacks

	V2: - Add max scale for multi-scale analysis
		- Exclude cholinergic bodies
*/

// Variables ---------------------------------------------------------------------
//	String variables
var dir								= "";										//	directory
var log_path						= "";										//	path for the log file
var micron							= getInfo("micrometer.abbreviation");		// 	micro symbol
var order							= "xyczt(default)";							//	hyperstack dimension order
var output_dir						= "";										//	dir for analysis output
var results							= "";										//	summarized results	

var spots_ch1_results				= "";										//	spot results ch 1
var spots_ch1_roi_set				= "";										//	spot ROI sets ch 1			
var spots_ch1_segmentation_method	= "Multi-Scale";								//	spots segmentation method ch 1
var spots_ch1_threshold				= "Fixed";									//	threshold method for spot segmentation ch 1
var spots_ch1_exclusion_threshold	= "Fixed";									//	threshold method for exclusion ch 1

var spots_ch2_results				= "";										//	spot results ch 2
var spots_ch2_roi_set				= "";										//	spot ROI sets ch 2	
var spots_ch2_segmentation_method	= "Multi-Scale";								//	spots segmentation method ch 2
var spots_ch2_threshold				= "Fixed";									//	threshold method for spot segmentation ch 2
var spots_ch2_exclusion_threshold	= "Fixed";									//	threshold method for exclusion ch 2

var spots_ch3_results				= "";										//	spot results ch 2
var spots_ch3_roi_set				= "";										//	spot ROI sets ch 3			
var spots_ch3_segmentation_method	= "Multi-Scale";								//	spots segmentation method ch 3
var spots_ch3_threshold				= "Fixed";									//	threshold method for spot segmentation ch 3
var spots_ch3_exclusion_threshold	= "Fixed";									//	threshold method for exclusion ch 3

var spots_ch4_results				= "";										//	spot results ch 3
var spots_ch4_roi_set				= "";										//	spot ROI sets ch 4	
var spots_ch4_segmentation_method	= "Multi-Scale";								//	spots segmentation method ch 4
var spots_ch4_threshold				= "Fixed";									//	threshold method for spot segmentation ch 4
var spots_ch4_exclusion_threshold	= "Fixed";									//	threshold method for exclusion ch 4

var suffix							= ".tif";									//	suffix for specifying the file type

//	Number variables
var channels						= 4;										//	number of channels
var image_height					= 1000;										//	image height
var image_width						= 1000;										//	image width
var pixel_size						= 0.1785703;								//	pixel size (µm)

var spots_ch1_filter_scale			= 2;										//	scale for Laplacian spots ch 1
var spots_ch1_fixed_threshold_value	= 1.5;										//	fixed maximum threshold for spot segmentation (if auto doesn't work well); ch 1
var spots_ch1_max_scale_multi		= 3; 										//  max scale multi-scale enhancement ch 1
var spots_ch1_min_area				= 0.2;										//	min spot size in micrometer ch 1
var spots_ch1_max_area				= 3;										//	max spot size in micrometer ch 1
var spots_ch1_exclusion_scale		= 15; 										//  scale medium filter to exclude cell bodies ch 1
var spots_ch1_exclusion_fixed_threshold	= 25; 									//  fixed threshold to exclude cell bodies ch 1
var spots_ch1_exclusion_min_area	= 20; 										//  min area of cell bodies to exclude in ch 1

var spots_ch2_filter_scale			= 2;										//	scale for Laplacian spots ch 2
var spots_ch2_fixed_threshold_value	= 3;										//	fixed maximum threshold for spot segmentation (if auto doesn't work well); ch 2
var spots_ch2_max_scale_multi		= 3; 										//  max scale multi-scale enhancement ch 2
var spots_ch2_max_area				= 3;										//	max spot size in micrometer ch 2
var spots_ch2_min_area				= 0.2;										//	min spot size in micrometer ch 2
var spots_ch2_exclusion_scale		= 15; 										//  scale medium filter to exclude cell bodies ch 2
var spots_ch2_exclusion_fixed_threshold	= 25; 									//  fixed threshold to exclude cell bodies ch 2
var spots_ch2_exclusion_min_area	= 20; 										//  min area of cell bodies to exclude in ch 2

var spots_ch3_filter_scale			= 1;										//	scale for Laplacian spots ch 3
var spots_ch3_fixed_threshold_value	= 2;										//	fixed maximum threshold for spot segmentation (if auto doesn't work well); ch 3
var spots_ch3_max_scale_multi		= 3; 										//  max scale multi-scale enhancement ch3
var spots_ch3_max_area				= 3;										//	max spot size in micrometer ch 3
var spots_ch3_min_area				= 0.2;										//	min spot size in micrometer ch 3
var spots_ch3_exclusion_scale		= 15; 										//  scale medium filter to exclude cell bodies ch 3
var spots_ch3_exclusion_fixed_threshold	= 25; 									//  fixed threshold to exclude cell bodies ch 3
var spots_ch3_exclusion_min_area	= 20; 										//  min area of cell bodies to exclude in ch 3

var spots_ch4_filter_scale			= 1;										//	scale for Laplacian spots ch 4
var spots_ch4_fixed_threshold_value	= 2;										//	fixed maximum threshold for spot segmentation (if auto doesn't work well); ch 4
var spots_ch4_max_scale_multi		= 3; 										//  max scale multi-scale enhancement ch4
var spots_ch4_max_area				= 3;										//	max spot size in micrometer ch 4
var spots_ch4_min_area				= 0.2;										//	min spot size in micrometer ch 4
var spots_ch4_exclusion_scale		= 15; 										//  scale medium filter to exclude cell bodies ch 4
var spots_ch4_exclusion_fixed_threshold	= 25; 									//  fixed threshold to exclude cell bodies ch 4
var spots_ch4_exclusion_min_area	= 20; 										//  min area of cell bodies to exclude in ch 4

var z_substack						= 3;
var z_step							= 2; 

//	Boolean Parameters
var spots_ch1						= false;
var spots_ch2						= true;
var spots_ch3						= true;
var spots_ch4						= true;

var spots_ch1_exclusion				= false;
var spots_ch2_exclusion				= false;
var spots_ch3_exclusion				= false;
var spots_ch4_exclusion				= true;

//	Arrays
var dimensions						= newArray("xyczt(default)","xyctz","xytcz","xytzc","xyztc","xyzct");		
var file_list						= newArray(0);								
var file_types 						= newArray(".tif",".tiff",".nd2",".ids",".jpg",".mvd2",".czi");		
var objects 						= newArray("spots_ch1","spots_ch2","spots_ch3","spots_ch4");
var prefixes 						= newArray(0);
var spot_segmentation_methods		= newArray("Gauss","Laplace","Multi-Scale");
var threshold_methods				= getList("threshold.methods");	
var threshold_methods				= Array.concat(threshold_methods,"Fixed");	

//---------------------------------------------------------------------Macros---------------------------------------------------------------------

macro "[I] Install Macro"{
	// Only works on my personal drive 
	run("Install...", "install=[/data/CBH/mverschuuren/ACAM/2021_BIL_MonicavandenBerg_Neurotransmitters/BIL_SynapseDetection_v2.ijm]");
}

macro "Autorun"{
	erase(1);
}

macro "Split Stack Into Fields Action Tool - C888 R0077 R9077 R9977 R0977"{
	setBatchMode(true);
	splitRegions();
	setBatchMode("exit and display");
}

macro "Z Project Action Tool - C888 R0099 R3399 R6699 R9999"{
	setBatchMode(true);
	projectImages();
	setBatchMode("exit and display");
}

macro "Setup Action Tool - C888 T5f16S"{
	setup();
}

macro "Analyse Single Image Action Tool - C888 T5f161"{
	setBatchMode(true);

	//Set Directory
	dir = getInfo("image.directory");
	output_dir = dir+"Output"+File.separator;
	if(!File.exists(output_dir)){
		print(output_dir);
		File.makeDirectory(output_dir);
	}
	start = getTime();

	//Image properties
	title = getTitle; 
	prefix = substring(title,0,lastIndexOf(title,suffix));
	setFileNames(prefix);
	id = getImageID;

	//Spot detection
	roiManager("reset");
	if(spots_ch1>0){
		args	= newArray(1,spots_ch1_segmentation_method,spots_ch1_threshold,spots_ch1_fixed_threshold_value,spots_ch1_filter_scale,spots_ch1_max_scale_multi,spots_ch1_min_area,spots_ch1_max_area,spots_ch1_exclusion,spots_ch1_exclusion_threshold,spots_ch1_exclusion_fixed_threshold,spots_ch1_exclusion_scale,spots_ch1_exclusion_min_area);
		snr 	= segmentSpots(id,args,false);
		if(snr>0){
			roiManager("Save",spots_ch1_roi_set);
			roiManager("reset");
		}
	}
	if(spots_ch2>0){
		args	= newArray(2,spots_ch2_segmentation_method,spots_ch2_threshold,spots_ch2_fixed_threshold_value,spots_ch2_filter_scale,spots_ch2_max_scale_multi,spots_ch2_min_area,spots_ch2_max_area,spots_ch2_exclusion,spots_ch2_exclusion_threshold,spots_ch2_exclusion_fixed_threshold,spots_ch2_exclusion_scale,spots_ch2_exclusion_min_area);
		snr 	= segmentSpots(id,args,false);
		if(snr>0)
		{
			roiManager("Save",spots_ch2_roi_set);
			roiManager("reset");
		}
	}
	if(spots_ch3>0){
		args	= newArray(3,spots_ch3_segmentation_method,spots_ch3_threshold,spots_ch3_fixed_threshold_value,spots_ch3_filter_scale,spots_ch3_max_scale_multi,spots_ch3_min_area,spots_ch3_max_area,spots_ch3_exclusion,spots_ch3_exclusion_threshold,spots_ch3_exclusion_fixed_threshold,spots_ch3_exclusion_scale,spots_ch3_exclusion_min_area);
		snr 	= segmentSpots(id,args,false);
		if(snr>0)
		{
			roiManager("Save",spots_ch3_roi_set);
			roiManager("reset");
		}
	}
	if(spots_ch4>0){
		args	= newArray(4,spots_ch4_segmentation_method,spots_ch4_threshold,spots_ch4_fixed_threshold_value,spots_ch4_filter_scale,spots_ch4_max_scale_multi,spots_ch4_min_area,spots_ch4_max_area,spots_ch4_exclusion,spots_ch4_exclusion_threshold,spots_ch4_exclusion_fixed_threshold,spots_ch4_exclusion_scale,spots_ch4_exclusion_min_area);
		snr 	= segmentSpots(id,args,false);
		if(snr>0)
		{
			roiManager("Save",spots_ch4_roi_set);
			roiManager("reset");
		}
	}

	//Spot analysis
	readout = analyzeRegions(id);
	if(readout){
		summaryResults();
	}
	
	print((getTime()-start)/1000,"sec");
	print("Analysis Done");
	setBatchMode("exit and display");
}

macro "Batch Analysis Action Tool - C888 T5f16#"{
	erase(1);
	setBatchMode(true);

	//Set Directory
	setDirectory();
	prefixes = scanFiles();
	fields = prefixes.length;
	setup();
	start = getTime();

	//Loop over images
	for(i=0;i<fields;i++){
		//Open image
		prefix = prefixes[i];
		file = prefix+suffix;
		setFileNames(prefix);
		print(i+1,"/",fields,":",prefix);
		path = dir+file;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		id 		= getImageID;

		//Spot detection
		roiManager("reset");
		if(spots_ch1>0){
			args	= newArray(1,spots_ch1_segmentation_method,spots_ch1_threshold,spots_ch1_fixed_threshold_value,spots_ch1_filter_scale,spots_ch1_max_scale_multi,spots_ch1_min_area,spots_ch1_max_area,spots_ch1_exclusion,spots_ch1_exclusion_threshold,spots_ch1_exclusion_fixed_threshold,spots_ch1_exclusion_scale,spots_ch1_exclusion_min_area);
			snr 	= segmentSpots(id,args,false);
			if(snr>0){
				roiManager("Save",spots_ch1_roi_set);
				roiManager("reset");
			}
		}
		if(spots_ch2>0){
			args	= newArray(2,spots_ch2_segmentation_method,spots_ch2_threshold,spots_ch2_fixed_threshold_value,spots_ch2_filter_scale,spots_ch2_max_scale_multi,spots_ch2_min_area,spots_ch2_max_area,spots_ch2_exclusion,spots_ch2_exclusion_threshold,spots_ch2_exclusion_fixed_threshold,spots_ch2_exclusion_scale,spots_ch2_exclusion_min_area);
			snr 	= segmentSpots(id,args,false);
			if(snr>0)
			{
				roiManager("Save",spots_ch2_roi_set);
				roiManager("reset");
			}
		}
		if(spots_ch3>0){
			args	= newArray(3,spots_ch3_segmentation_method,spots_ch3_threshold,spots_ch3_fixed_threshold_value,spots_ch3_filter_scale,spots_ch3_max_scale_multi,spots_ch3_min_area,spots_ch3_max_area,spots_ch3_exclusion,spots_ch3_exclusion_threshold,spots_ch3_exclusion_fixed_threshold,spots_ch3_exclusion_scale,spots_ch3_exclusion_min_area);
			snr 	= segmentSpots(id,args,false);
			if(snr>0)
			{
				roiManager("Save",spots_ch3_roi_set);
				roiManager("reset");
			}
		}
		if(spots_ch4>0){
			args	= newArray(4,spots_ch4_segmentation_method,spots_ch4_threshold,spots_ch4_fixed_threshold_value,spots_ch4_filter_scale,spots_ch4_max_scale_multi,spots_ch4_min_area,spots_ch4_max_area,spots_ch4_exclusion,spots_ch4_exclusion_threshold,spots_ch4_exclusion_fixed_threshold,spots_ch4_exclusion_scale,spots_ch4_exclusion_min_area);
			snr 	= segmentSpots(id,args,false);
			if(snr>0)
			{
				roiManager("Save",spots_ch4_roi_set);
				roiManager("reset");
			}
		}

		//Spot analysis
		readout = analyzeRegions(id);
		selectImage(id); close;
		if(readout){
			summaryResults();
		}
		erase(0);
	}

	//Concatenate results
	concatenateResults();
	saveAs(".txt",output_dir+"ConcatenatedResults.txt");
	
	print((getTime()-start)/1000,"sec");
	if(isOpen("Log")){
		selectWindow("Log");
		saveAs("txt",log_path);
	}
	print("Complete Analysis Done");
	setBatchMode("exit and display");
}

macro "Segment Spots Action Tool - C888 H00f5f8cf3f0800 Cf88 o3222 o4b22 o5822 o7522 oa822"
{
	erase(0);
	setBatchMode(true);
	id 		= getImageID;
	c 		= getNumber("Spot Channel",1);
	if(c == 1){
		args = newArray(1,spots_ch1_segmentation_method,spots_ch1_threshold,spots_ch1_fixed_threshold_value,spots_ch1_filter_scale,spots_ch1_max_scale_multi,spots_ch1_min_area,spots_ch1_max_area,spots_ch1_exclusion,spots_ch1_exclusion_threshold,spots_ch1_exclusion_fixed_threshold,spots_ch1_exclusion_scale,spots_ch1_exclusion_min_area);
	}
	if(c == 2){
		args = newArray(2,spots_ch2_segmentation_method,spots_ch2_threshold,spots_ch2_fixed_threshold_value,spots_ch2_filter_scale,spots_ch2_max_scale_multi,spots_ch2_min_area,spots_ch2_max_area,spots_ch2_exclusion,spots_ch2_exclusion_threshold,spots_ch2_exclusion_fixed_threshold,spots_ch2_exclusion_scale,spots_ch2_exclusion_min_area);
	}
	if(c == 3){
		args = newArray(3,spots_ch3_segmentation_method,spots_ch3_threshold,spots_ch3_fixed_threshold_value,spots_ch3_filter_scale,spots_ch3_max_scale_multi,spots_ch3_min_area,spots_ch3_max_area,spots_ch3_exclusion,spots_ch3_exclusion_threshold,spots_ch3_exclusion_fixed_threshold,spots_ch3_exclusion_scale,spots_ch3_exclusion_min_area);
	}
	if(c == 4){
		args = newArray(4,spots_ch4_segmentation_method,spots_ch4_threshold,spots_ch4_fixed_threshold_value,spots_ch4_filter_scale,spots_ch4_max_scale_multi,spots_ch4_min_area,spots_ch4_max_area,spots_ch4_exclusion,spots_ch4_exclusion_threshold,spots_ch4_exclusion_fixed_threshold,spots_ch4_exclusion_scale,spots_ch4_exclusion_min_area);
	}
	snr 	= segmentSpots(id,args,true);
	selectImage(id);
	toggleOverlay();
	setBatchMode("exit and display");
	run("Tile");
	

}

macro "Toggle Overlay Action Tool - Caaa O11ee"{
	toggleOverlay();
}

macro "[t] Toggle Overlay"{
	toggleOverlay();
}

macro "Verification Stack Action Tool - C888 T5f16V"{
	erase(1);
	setBatchMode(true);
	setDirectory();
	prefixes = scanFiles();
	names = prefixes;
	createOverlay(names);
	setBatchMode("exit and display");
	run("Channels Tool... ");
}

macro "Concatenate Results Action Tool - C888 T5f16C"
{
	setBatchMode(true);
	setDirectory();
	concatenateResults();
	saveAs(".txt",output_dir+"ConcatenatedResults.txt");
	setBatchMode("exit and display");
}

//---------------------------------------------------------------------Functions---------------------------------------------------------------------

function splitRegions(){
	erase(1);

	//GUI
	Dialog.create("Split Fields...");
	Dialog.addString("Destination Directory Name","Export",25);
	Dialog.addString("Add a prefix","",25);
	Dialog.addChoice("Import format",file_types,".mvd2");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.addNumber("Channels",4);
	Dialog.addChoice("Dimension order",dimensions, order);
	Dialog.show;
	dest 		= Dialog.getString;
	pre			= Dialog.getString;
	ext			= Dialog.getChoice;
	suffix 		= Dialog.getChoice;
	channels	= Dialog.getNumber;
	order 		= Dialog.getChoice;

	//Set directory
	dir = getDirectory("");
	file_list = getFileList(dir);
	destination_dir = dir+dest+File.separator;
	File.makeDirectory(destination_dir);

	//Loop over images: open and save individual files
	for(i=0;i<file_list.length;i++){
		path = dir+file_list[i];
		if(endsWith(path,ext)){		
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_all_series view=Hyperstack ");
			while(nImages>0){
				selectImage(nImages);
				id=getImageID();
				selectImage(id);
				run("Stack to Hyperstack...", "order="+order+" channels="+channels+" slices="+(nSlices/channels)+" frames=1 display=Color");
				id=getImageID();
				title = getTitle; 
				saveAs(suffix,destination_dir+title+suffix);
				selectImage(id); close;
			}
		}
	}
	print("Done");
}


function projectImages(){
	erase(1);

	//GUI
	Dialog.create("Project Images...");
	Dialog.addChoice("Import format",file_types,".tif");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.addNumber("Z step", z_step, 0, 2, micron);
	Dialog.addNumber("# Slices of Substack (0=all)", z_substack, 0, 2, "");
	Dialog.show;
	ext			= Dialog.getChoice;
	suffix 		= Dialog.getChoice;
	z_step		= Dialog.getNumber();
	z_substack	= Dialog.getNumber();
	
	//Set directory
	dir 		= getDirectory("");
	file_list 	= getFileList(dir);
	folder_subSample = "SubSample";
	dest_subSample 	= dir+folder_subSample+File.separator;
	File.makeDirectory(dest_subSample);
	folder_max = "MaxProj";
	dest_max	= dest_subSample+folder_max+File.separator;
	File.makeDirectory(dest_max);

	//Loop over images: substack if needed+ max projection
	for(i=0;i<file_list.length;i++){
		path = dir+file_list[i];
		if(endsWith(path,ext)){		
			print(i+1);

			//Open image + check properties
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default concatenate_series open_all_series view=Hyperstack ");
			ns = nSlices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+ns/channels+" frames=1 display=Color");
			id = getImageID;
			title = getTitle;
			print(title);
			getDimensions(width, height, channels, slices, frames);
			getVoxelSize(width, height, depth, unit);
			print("Res. z: " + depth + " " + unit);

			//Subsample if z-resolution < z_step
			run("Duplicate...", "duplicate");
			idDup=getImageID();
			if(unit!="microns"){
				showMessageWithCancel("Maximum intensity projection","Add image calibration in microns");
				selectImage(id); close;
				selectImage(idDup); close;
				break;
			}else{
				if(z_step>2){
					showMessageWithCancel("Maximum intensity projection","Z-range not included in script");
					selectImage(id); close;
					selectImage(idDup); close;
					break;
				}
				if(depth>z_step){
					showMessageWithCancel("Maximum intensity projection","Voxel depth larger than z_step");
					selectImage(id); close;
					selectImage(idDup); close;
					break;
				}
				if(depth<z_step){
					factor=z_step/depth;
					s=2; 
					selectImage(idDup);
					getDimensions(width, height, channels, slices, frames);
					finalSlices=Math.ceil(slices/factor);
					print("Number of slices before subsampling: "+ slices);
					while(slices>finalSlices){
						selectImage(idDup);
						Stack.setSlice(s);
						run("Delete Slice", "delete=slice");
						selectImage(idDup);
						getDimensions(width, height, channels, slices, frames);
						s=s+1;
					}
					selectImage(idDup);
					getVoxelSize(width, height, depth, unit);
					setVoxelSize(width, height, z_step, unit);
					selectImage(idDup);
					getDimensions(width, height, channels, slices, frames);
					print("Number of slices after subsampling: "+ slices);
				}
			}

			selectImage(idDup); 
			saveAs(suffix,dest_subSample+"SubSample_"+title+suffix);

			//Addition substack if needed
			if(z_substack!=0){
				selectImage(idDup);
				getDimensions(width, height, channels, slices, frames);
				print("Slices in stack: " + slices );
				if(slices>z_substack){
					zMin=round((slices-z_substack)/2)+1;
					zMax=zMin+z_substack-1;
					selectImage(idDup);
					run("Duplicate...", "duplicate slices="+zMin+"-"+zMax);
					idDupSub=getImageID();
					getDimensions(width, height, channels, slices, frames);
					print(slices +" slices in substack: "+zMin+" - "+zMax);
					run("Z Project...", "projection=[Max Intensity]");
					zid = getImageID;	
					selectImage(zid); saveAs(suffix,dest_max+"Max_"+title+suffix);
					selectImage(idDupSub); close;	
					selectImage(zid); close;
				}else{
					run("Z Project...", "projection=[Max Intensity]");
					zid = getImageID;
					selectImage(zid); saveAs(suffix,dest_max+"Max_"+title+suffix);
					selectImage(zid); close;	
				}
			}else{
				run("Z Project...", "projection=[Max Intensity]");
				zid = getImageID;
				selectImage(zid); saveAs(suffix,dest_max+"Max_"+title+suffix);
				selectImage(zid); close;	
			}
			selectImage(id); close;
			selectImage(idDup); close;		
		}
	}
	print("Done");
}

function setOptions(){
	run("Options...", "iterations=1 count=1");
	run("Colors...", "foreground=white correct_background=black selection=yellow");
	run("Overlay Options...", "stroke=red width=1 fill=none");
	setBackgroundColor(0, 0, 0);
	setForegroundColor(255,255,255);
}

function getMoment(){
     MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+second;
     return TimeString;
}

function erase(all){
	if(all){print("\\Clear");run("Close All");}
	run("Clear Results");
	roiManager("reset");
	run("Collect Garbage");
}

function setDirectory(){
	dir = getDirectory("Choose a Source Directory");
	file_list = getFileList(dir);
	output_dir = dir+"Output"+File.separator;
	if(!File.exists(output_dir))File.makeDirectory(output_dir);
	log_path = output_dir+"Log.txt";
}

function setFileNames(prefix){
	spots_ch1_roi_set	= output_dir+prefix+"_spots_ch1_roi_set.zip";
	spots_ch2_roi_set	= output_dir+prefix+"_spots_ch2_roi_set.zip";
	spots_ch3_roi_set	= output_dir+prefix+"_spots_ch3_roi_set.zip";
	spots_ch4_roi_set	= output_dir+prefix+"_spots_ch4_roi_set.zip";
	spots_ch1_results	= output_dir+prefix+"_spots_ch1_results.txt";
	spots_ch2_results	= output_dir+prefix+"_spots_ch2_results.txt";
	spots_ch3_results	= output_dir+prefix+"_spots_ch3_results.txt";
	spots_ch4_results	= output_dir+prefix+"_spots_ch4_results.txt";
	results				= output_dir+prefix+"_summary.txt";
}

function scanFiles(){
	prefixes = newArray(0);
	for(i=0;i<file_list.length;i++){
		path = dir+file_list[i];
		if(endsWith(path,suffix) && indexOf(path,"flatfield")<0){
			print(path);
			prefixes = Array.concat(prefixes,substring(file_list[i],0,lastIndexOf(file_list[i],suffix)));			
		}
	}
	return prefixes;
}

function setup(){
	setOptions();
	Dialog.create("Colocalisation neurotransmitters");
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Image Type", file_types, suffix);
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Pixel Size", pixel_size, 3, 5, micron+" (only if not calibrated)");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Number of Channels", channels, 0, 5, " ");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Spot channels:");
	Dialog.setInsets(0,0,0);
	labels = newArray(4);defaults = newArray(4);
	labels[0] = "Channel 1";		defaults[0] = spots_ch1;
	labels[1] = "Channel 2";		defaults[1] = spots_ch2;
	labels[2] = "Channel 3";		defaults[2] = spots_ch3;
	labels[3] = "Channel 4";		defaults[3] = spots_ch4;
	Dialog.setInsets(0,0,0);
	Dialog.addCheckboxGroup(1,4,labels,defaults);
	Dialog.setInsets(0,0,0);
	Dialog.show();
	
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
	
	suffix							= Dialog.getChoice();		print("Image type:",suffix);
	pixel_size						= Dialog.getNumber(); 		print("Pixel size:",pixel_size);
	channels 						= Dialog.getNumber();		print("Channels:",channels);
	spots_ch1						= Dialog.getCheckbox(); 	print("Spots channel 1",spots_ch1);
	spots_ch2	 					= Dialog.getCheckbox();		print("Spots channel 2",spots_ch2);
	spots_ch3	 					= Dialog.getCheckbox();		print("Spots channel 3:",spots_ch3);
	spots_ch4						= Dialog.getCheckbox();		print("Spots channel 4:",spots_ch4);

	Dialog.create("Spot Analysis");
	if(spots_ch1){
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------------------------------------------------------------------------------  Channel 1  -----------------------------------------------------------------------------------------------", 14, "#ff0000");
		Dialog.setInsets(0,0,0);	
		Dialog.addChoice("Spot Segmentation Method", spot_segmentation_methods, spots_ch1_segmentation_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Gauss/Laplace Scale", spots_ch1_filter_scale, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Max Scale Multi", spots_ch1_max_scale_multi, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Threshold Algorithm", threshold_methods, spots_ch1_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", spots_ch1_fixed_threshold_value, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Minimum Spot Size", spots_ch1_min_area, 2, 4, micron);
		Dialog.addToSameRow();
		Dialog.addNumber("Maximum Spot Size", spots_ch1_max_area, 2, 4, micron);
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("");
		Dialog.setInsets(0,0,0);
		Dialog.addCheckbox("Exclusion cell bodies", spots_ch1_exclusion);
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Exclusion: Threshold Algorithm", threshold_methods, spots_ch1_exclusion_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Fixed Threshold", spots_ch1_exclusion_fixed_threshold, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Median filter Scale", spots_ch1_exclusion_scale, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Minimum area", spots_ch1_exclusion_min_area, 2, 4, micron);
	}
	if(spots_ch2){
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------------------------------------------------------------------------------  Channel 2  -----------------------------------------------------------------------------------------------", 14, "#ff0000");
		Dialog.setInsets(0,0,0);	
		Dialog.addChoice("Spot Segmentation Method", spot_segmentation_methods, spots_ch2_segmentation_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Gauss/Laplace Scale", spots_ch2_filter_scale, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Max Scale Multi", spots_ch2_max_scale_multi, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Threshold Algorithm", threshold_methods, spots_ch2_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", spots_ch2_fixed_threshold_value, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Minimum Spot Size", spots_ch2_min_area, 2, 4, micron);
		Dialog.addToSameRow();
		Dialog.addNumber("Maximum Spot Size", spots_ch2_max_area, 2, 4, micron);
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("");
		Dialog.setInsets(0,0,0);
		Dialog.addCheckbox("Exclusion cell bodies", spots_ch2_exclusion);
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Exclusion: Threshold Algorithm", threshold_methods, spots_ch2_exclusion_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Fixed Threshold", spots_ch2_exclusion_fixed_threshold, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Median filter Scale", spots_ch2_exclusion_scale, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Minimum area", spots_ch2_exclusion_min_area, 2, 4, micron);
	}
	if(spots_ch3){
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------------------------------------------------------------------------------  Channel 3  -----------------------------------------------------------------------------------------------", 14, "#ff0000");
		Dialog.setInsets(0,0,0);	
		Dialog.addChoice("Spot Segmentation Method", spot_segmentation_methods, spots_ch3_segmentation_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Gauss/Laplace Scale", spots_ch3_filter_scale, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Max Scale Multi", spots_ch3_max_scale_multi, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Threshold Algorithm", threshold_methods, spots_ch3_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", spots_ch3_fixed_threshold_value, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Minimum Spot Size", spots_ch3_min_area, 2, 4, micron);
		Dialog.addToSameRow();
		Dialog.addNumber("Maximum Spot Size", spots_ch3_max_area, 2, 4, micron);
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("");
		Dialog.setInsets(0,0,0);
		Dialog.addCheckbox("Exclusion cell bodies", spots_ch3_exclusion);
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Exclusion: Threshold Algorithm", threshold_methods, spots_ch3_exclusion_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Fixed Threshold", spots_ch3_exclusion_fixed_threshold, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Median filter Scale", spots_ch3_exclusion_scale, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Minimum area", spots_ch3_exclusion_min_area, 2, 4, micron);
	}
	if(spots_ch4){
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------------------------------------------------------------------------------  Channel 4  -----------------------------------------------------------------------------------------------", 14, "#ff0000");
		Dialog.setInsets(0,0,0);	
		Dialog.addChoice("Spot Segmentation Method", spot_segmentation_methods, spots_ch4_segmentation_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Gauss/Laplace Scale", spots_ch4_filter_scale, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Max Scale Multi", spots_ch4_max_scale_multi, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Threshold Algorithm", threshold_methods, spots_ch4_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", spots_ch4_fixed_threshold_value, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Minimum Spot Size", spots_ch4_min_area, 2, 4, micron);
		Dialog.addToSameRow();
		Dialog.addNumber("Maximum Spot Size", spots_ch4_max_area, 2, 4, micron);
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("");
		Dialog.setInsets(0,0,0);
		Dialog.addCheckbox("Exclusion cell bodies", spots_ch4_exclusion);
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Exclusion: Threshold Algorithm", threshold_methods, spots_ch4_exclusion_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Fixed Threshold", spots_ch4_exclusion_fixed_threshold, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Median filter Scale", spots_ch4_exclusion_scale, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Exclusion: Minimum area", spots_ch4_exclusion_min_area, 2, 4, micron);
	}
	Dialog.show();

	if(spots_ch1){
		spots_ch1_segmentation_method		= Dialog.getChoice();		print("CH1 - Spot Segmentation Method:",spots_ch1_segmentation_method);
		spots_ch1_filter_scale 				= Dialog.getNumber();		print("CH1 - Laplace scale:",spots_ch1_filter_scale);
		spots_ch1_max_scale_multi			= Dialog.getNumber();		print("CH1 - Max multi scale:",spots_ch1_max_scale_multi);
		spots_ch1_threshold 				= Dialog.getChoice();		print("CH1 - Spot autoThreshold:",spots_ch1_threshold);
		spots_ch1_fixed_threshold_value  	= Dialog.getNumber();		print("CH1 - Fixed threshold:",spots_ch1_fixed_threshold_value);
		spots_ch1_min_area	 				= Dialog.getNumber();		print("CH1 - Min. spot size:", spots_ch1_min_area);
		spots_ch1_max_area 					= Dialog.getNumber();		print("CH1 - Max. spot size:",spots_ch1_max_area);
		spots_ch1_exclusion					= Dialog.getCheckbox();		print("CH1 - Exclusion cell bodies:",spots_ch1_exclusion);
		spots_ch1_exclusion_threshold		= Dialog.getChoice();		print("CH1 - Exclusion auto threshold:",spots_ch1_exclusion_threshold);
		spots_ch1_exclusion_fixed_threshold	= Dialog.getNumber();		print("CH1 - Exclusion fixed threshold:",spots_ch1_exclusion_fixed_threshold);
		spots_ch1_exclusion_scale			= Dialog.getNumber();		print("CH1 - Exclusion median filter scale:", spots_ch1_exclusion_scale); 
		spots_ch1_exclusion_min_area		= Dialog.getNumber();		print("CH1 - Exclusion min area:", spots_ch1_exclusion_min_area);
	}
	if(spots_ch2){
		spots_ch2_segmentation_method		= Dialog.getChoice();		print("CH2 - Spot Segmentation Method:",spots_ch2_segmentation_method);
		spots_ch2_filter_scale 				= Dialog.getNumber();		print("CH2 - Laplace scale:",spots_ch2_filter_scale);
		spots_ch2_max_scale_multi			= Dialog.getNumber();		print("CH2 - Max multi scale:",spots_ch2_max_scale_multi);
		spots_ch2_threshold 				= Dialog.getChoice();		print("CH2 - Spot autoThreshold:",spots_ch2_threshold);
		spots_ch2_fixed_threshold_value  	= Dialog.getNumber();		print("CH2 - Fixed threshold:",spots_ch2_fixed_threshold_value);
		spots_ch2_min_area	 				= Dialog.getNumber();		print("CH2 - Min. spot size:", spots_ch2_min_area);
		spots_ch2_max_area 					= Dialog.getNumber();		print("CH2 - Max. spot size:",spots_ch2_max_area);
		spots_ch2_exclusion					= Dialog.getCheckbox();		print("CH2 - Exclusion cell bodies:",spots_ch2_exclusion);
		spots_ch2_exclusion_threshold		= Dialog.getChoice();		print("CH2 - Exclusion auto threshold:",spots_ch2_exclusion_threshold);
		spots_ch2_exclusion_fixed_threshold	= Dialog.getNumber();		print("CH2 - Exclusion fixed threshold:",spots_ch2_exclusion_fixed_threshold);
		spots_ch2_exclusion_scale			= Dialog.getNumber();		print("CH2 - Exclusion median filter scale:", spots_ch2_exclusion_scale); 
		spots_ch2_exclusion_min_area		= Dialog.getNumber();		print("CH2 - Exclusion min area:", spots_ch2_exclusion_min_area);
	}
	if(spots_ch3){
		spots_ch3_segmentation_method		= Dialog.getChoice();		print("CH3 - Spot Segmentation Method:",spots_ch3_segmentation_method);
		spots_ch3_filter_scale 				= Dialog.getNumber();		print("CH3 - Laplace scale:",spots_ch3_filter_scale);
		spots_ch3_max_scale_multi			= Dialog.getNumber();		print("CH3 - Max multi scale:",spots_ch3_max_scale_multi);
		spots_ch3_threshold 				= Dialog.getChoice();		print("CH3 - Spot autoThreshold:",spots_ch3_threshold);
		spots_ch3_fixed_threshold_value  	= Dialog.getNumber();		print("CH3 - Fixed threshold:",spots_ch3_fixed_threshold_value);
		spots_ch3_min_area	 				= Dialog.getNumber();		print("CH3 - Min. spot size:", spots_ch3_min_area);
		spots_ch3_max_area 					= Dialog.getNumber();		print("CH3 - Max. spot size:",spots_ch3_max_area);
		spots_ch3_exclusion					= Dialog.getCheckbox();		print("CH3 - Exclusion cell bodies:",spots_ch3_exclusion);
		spots_ch3_exclusion_threshold		= Dialog.getChoice();		print("CH3 - Exclusion auto threshold:",spots_ch3_exclusion_threshold);
		spots_ch3_exclusion_fixed_threshold	= Dialog.getNumber();		print("CH3 - Exclusion fixed threshold:",spots_ch3_exclusion_fixed_threshold);
		spots_ch3_exclusion_scale			= Dialog.getNumber();		print("CH3 - Exclusion median filter scale:", spots_ch3_exclusion_scale); 
		spots_ch3_exclusion_min_area		= Dialog.getNumber();		print("CH3 - Exclusion min area:", spots_ch3_exclusion_min_area);
	}
	if(spots_ch4){
		spots_ch4_segmentation_method		= Dialog.getChoice();		print("CH4 - Spot Segmentation Method:",spots_ch4_segmentation_method);
		spots_ch4_filter_scale 				= Dialog.getNumber();		print("CH4 - Laplace scale:",spots_ch4_filter_scale);
		spots_ch4_max_scale_multi			= Dialog.getNumber();		print("CH4 - Max multi scale:",spots_ch4_max_scale_multi);		
		spots_ch4_threshold 				= Dialog.getChoice();		print("CH4 - Spot autoThreshold:",spots_ch4_threshold);
		spots_ch4_fixed_threshold_value  	= Dialog.getNumber();		print("CH4 - Fixed threshold:",spots_ch4_fixed_threshold_value);
		spots_ch4_min_area	 				= Dialog.getNumber();		print("CH4 - Min. spot size:", spots_ch4_min_area);
		spots_ch4_max_area 					= Dialog.getNumber();		print("CH4 - Max. spot size:",spots_ch4_max_area);
		spots_ch4_exclusion					= Dialog.getCheckbox();		print("CH4 - Exclusion cell bodies:",spots_ch4_exclusion);
		spots_ch4_exclusion_threshold		= Dialog.getChoice();		print("CH4 - Exclusion auto threshold:",spots_ch4_exclusion_threshold);
		spots_ch4_exclusion_fixed_threshold	= Dialog.getNumber();		print("CH4 - Exclusion fixed threshold:",spots_ch4_exclusion_fixed_threshold);
		spots_ch4_exclusion_scale			= Dialog.getNumber();		print("CH4 - Exclusion median filter scale:", spots_ch4_exclusion_scale); 
		spots_ch4_exclusion_min_area		= Dialog.getNumber();		print("CH4 - Exclusion min area:", spots_ch4_exclusion_min_area);	
	}
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
}

function calibrateImage(id)
{
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(unit!=micron)run("Properties...", " unit="+micron+" pixel_width="+pixel_size+" pixel_height="+pixel_size);
	else pixel_size = pixelWidth;
}

function decalibrateImage(id)
{
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(unit!="pixel")run("Properties...", " unit=pixel pixel_width=1 pixel_height=1");
}

function segmentSpots(id,args,booleanDebug){
	//Extract features
	spot_channel				= args[0];
	spot_segmentation_method	= args[1];
	spot_threshold_method		= args[2];
	spot_fixed_threshold_value	= args[3];
	scale						= args[4];
	max_scale_multi				= args[5];
	spot_min_area				= args[6];
	spot_max_area				= args[7];
	exclusion 					= args[8];	
	exclusion_segmentation_method = args[9];
	exclusion_fixed_threshold 	= args[10];
	exclusion_scale 			= args[11];
	exclusion_min_area 			= args[12];

	//Image properties
	selectImage(id);
	run("Select None");
	if(Stack.isHyperstack){
		run("Duplicate...", "title=copy duplicate channels="+spot_channel);	
	}
	else{
		setSlice(spot_channel);
		run("Duplicate...","title=copy ");
	}
	cid = getImageID;
	decalibrateImage(cid);
	selectImage(cid);
	title = getTitle;

	//Preprocessing
	if(spot_segmentation_method=="Gauss"){
		run("Duplicate...","title=Gauss");
		run("Gaussian Blur...", "sigma="+scale);
		run("Invert");
		idLap = getImageID;
	}
	else if(spot_segmentation_method=="Laplace"){
		run("FeatureJ Laplacian", "compute smoothing="+scale);
		idLap = getImageID;
	}
	else if(spot_segmentation_method=="Multi-Scale") {
		cid=getImageID;
		title=getTitle;
		scale=max_scale_multi;
		e=0;
		while(e<scale){			
			e++;
			selectImage(cid);
			run("FeatureJ Laplacian", "compute smoothing="+e);
			selectWindow(title+" Laplacian");
			run("Multiply...","value="+e*e);
			rename("scale "+e);
			eid = getImageID;
			if(e>1)
			{
				selectImage(eid);run("Select All");run("Copy");close;
				selectWindow("scale 1");run("Add Slice");run("Paste");
			}
		}
		selectWindow("scale 1"); 
		nlid = getImageID;
		selectImage(nlid);
		run("Z Project...", "start=1 projection=[Sum Slices]");
		idLap = getImageID;
		selectImage(nlid); close;
	}

	//Thresholding
	selectImage(idLap);
	run("Duplicate...","title=lapBinary");
	idLapBinary=getImageID();	
	selectImage(idLapBinary);
	setOption("BlackBackground",false);
	if(spot_threshold_method=="Fixed"){
		setAutoThreshold("Default ");
		getThreshold(mit,mat); 
		print("Original threshold settings:",mit,mat);
		setThreshold(minOf(mit,-spot_fixed_threshold_value),-spot_fixed_threshold_value);
		print("Fixed threshold settings:",minOf(mit,-spot_fixed_threshold_value),-spot_fixed_threshold_value);
	}
	else{
		setAutoThreshold(spot_threshold_method+" ");
	}
	getThreshold(mit,mat); 
	print("Threshold:",mit,mat);
	run("Convert to Mask");

	//Max finding
	selectImage(idLap);
	diameter=sqrt(spot_min_area/PI)*2;
	run("Find Maxima...", "prominence="+diameter+" light output=[Segmented Particles]");
	idRegions=getImageID();
	rename("Regions");
	imageCalculator("AND","lapBinary", "Regions");

	//Exclude cell bodies
	if(exclusion){
		selectImage(cid);
		run("Duplicate...", " ");
		idExcl=getImageID();
		run("Median...", "radius="+exclusion_scale);
		
		if(exclusion_segmentation_method=="Fixed"){
			getStatistics(area, mean, min, max, std, histogram);
			setAutoThreshold("Default ");
			setThreshold(exclusion_fixed_threshold,max);
		}
		else{
			setAutoThreshold(exclusion_segmentation_method+" ");
		}
		run("Convert to Mask");	
		calibrateImage(idExcl);
		run("Analyze Particles...", "size="+exclusion_min_area+"-Infinity show=Masks");
		idMaskExcl=getImageID();
		rename("Exclude-CH"+spot_channel);
		run("Invert");
		imageCalculator("AND","lapBinary", "Exclude-CH"+spot_channel);
		selectImage(idExcl); close();
		if(!booleanDebug){
			selectImage(idMaskExcl); close();
		}
	}

	//Measurements
	selectImage(idLapBinary);
	calibrateImage(idLapBinary);
	run("Set Measurements...", "  area min mean redirect=None decimal=4");
	run("Analyze Particles...", "size="+spot_min_area+"-"+spot_max_area+" circularity=0.00-1.00 show=Nothing display clear include add");
	snr = roiManager("count"); 
	print(snr,"spots");
	
	//To avoid excessive spot finding when there are no true spots
	if(snr>10000){
		print("excessive number, hence reset"); 
		snr=0; 
		roiManager("reset");
	}
	selectImage(idLap); close;
	selectImage(cid); close;
	selectImage(idLapBinary); close;
	selectImage(idRegions); close;
	
	return snr;
}

function analyzeRegions(id)
{
	erase(0); 
	mask = 0;
	readout = 1;
	selectImage(id);
	calibrateImage(id);
	//	analyze spot rois
	if(File.exists(spots_ch1_roi_set)){
		ms = nSlices;
		run("Set Measurements...", "  area mean min redirect=None decimal=4");
		roiManager("Open",spots_ch1_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		run("Clear Results");
		IJ.renameResults("Temp","Results");
		updateResults;
		saveAs("Measurements",spots_ch1_results);
		erase(0);
	}
	if(File.exists(spots_ch2_roi_set)){
		ms = nSlices;
		run("Set Measurements...", "  area mean min redirect=None decimal=4");
		roiManager("Open",spots_ch2_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		run("Clear Results");
		IJ.renameResults("Temp","Results");
		updateResults;
		saveAs("Measurements",spots_ch2_results);
		erase(0);
	}
	if(File.exists(spots_ch3_roi_set)){
		ms = nSlices;
		run("Set Measurements...", "  area mean min redirect=None decimal=4");
		roiManager("Open",spots_ch3_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		run("Clear Results");
		IJ.renameResults("Temp","Results");
		updateResults;
		saveAs("Measurements",spots_ch3_results);
		erase(0);
	}
	if(File.exists(spots_ch4_roi_set)){
		ms = nSlices;
		run("Set Measurements...", "  area mean min redirect=None decimal=4");
		roiManager("Open",spots_ch4_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		run("Clear Results");
		IJ.renameResults("Temp","Results");
		updateResults;
		saveAs("Measurements",spots_ch4_results);
		erase(0);
	}
	return readout;
}

function summaryResults(){
	erase(0); 
	row=0;
	resultWindow=false;
	
	if(File.exists(spots_ch1_results)){
		if(isOpen("Results")){
			resultWindow=true;
			IJ.renameResults("Results","Temp");
		}
		run("Results... ", "open=["+spots_ch1_results+"]");
		count = nResults;
		resultLabels 	= getResultLabels();
		matrix 			= results2matrix(resultLabels);
		if(resultWindow){
			selectWindow("Results"); run("Close");
			IJ.renameResults("Temp","Results");
		}else{
			run("Clear Results");
		}
		for(s=0;s<resultLabels.length;s++){
			value=0;
			for(r=0;r<count;r++){
				selectImage(matrix);
				p = getPixel(s,r);
				value+=p;
			}
			setResult("Spot_SC"+1+"_Nr",0,count);
			setResult("Spot_SC"+1+"_"+resultLabels[s]+"_Sum",0,value);              
			setResult("Spot_SC"+1+"_"+resultLabels[s]+"_Mean",0,value/count);
					
		}
		selectImage(matrix); close;
		updateResults();
		resultWindow=false;
	}

	if(File.exists(spots_ch2_results)){
		if(isOpen("Results")){
			resultWindow=true;
			IJ.renameResults("Results","Temp");
		}
		run("Results... ", "open=["+spots_ch2_results+"]");
		count = nResults;
		resultLabels 	= getResultLabels();
		matrix 			= results2matrix(resultLabels);
		if(resultWindow){
			selectWindow("Results"); run("Close");
			IJ.renameResults("Temp","Results");
		}else{
			run("Clear Results");
		}
		for(s=0;s<resultLabels.length;s++){
			value=0;
			for(r=0;r<count;r++){
				selectImage(matrix);
				p = getPixel(s,r);
				value+=p;
			}
			setResult("Spot_SC"+2+"_Nr",0,count);
			setResult("Spot_SC"+2+"_"+resultLabels[s]+"_Sum",0,value);              
			setResult("Spot_SC"+2+"_"+resultLabels[s]+"_Mean",0,value/count);
					
		}
		selectImage(matrix); close;
		updateResults();
		resultWindow=false;
	}

	if(File.exists(spots_ch3_results)){
		if(isOpen("Results")){
			resultWindow=true;
			IJ.renameResults("Results","Temp");
		}
		run("Results... ", "open=["+spots_ch3_results+"]");
		count = nResults;
		resultLabels 	= getResultLabels();
		matrix 			= results2matrix(resultLabels);
		if(resultWindow){
			selectWindow("Results"); run("Close");
			IJ.renameResults("Temp","Results");
		}else{
			run("Clear Results");
		}
		for(s=0;s<resultLabels.length;s++){
			value=0;
			for(r=0;r<count;r++){
				selectImage(matrix);
				p = getPixel(s,r);
				value+=p;
			}
			setResult("Spot_SC"+3+"_Nr",0,count);
			setResult("Spot_SC"+3+"_"+resultLabels[s]+"_Sum",0,value);              
			setResult("Spot_SC"+3+"_"+resultLabels[s]+"_Mean",0,value/count);
					
		}
		selectImage(matrix); close;
		updateResults();
		resultWindow=false;
	}

	if(File.exists(spots_ch4_results)){
		if(isOpen("Results")){
			resultWindow=true;
			IJ.renameResults("Results","Temp");
		}
		run("Results... ", "open=["+spots_ch4_results+"]");
		count = nResults;
		resultLabels 	= getResultLabels();
		matrix 			= results2matrix(resultLabels);
		if(resultWindow){
			selectWindow("Results"); run("Close");
			IJ.renameResults("Temp","Results");
		}else{
			run("Clear Results");
		}
		for(s=0;s<resultLabels.length;s++){
			value=0;
			for(r=0;r<count;r++){
				selectImage(matrix);
				p = getPixel(s,r);
				value+=p;
			}
			setResult("Spot_SC"+4+"_Nr",row,count);
			setResult("Spot_SC"+4+"_"+resultLabels[s]+"_Sum",0,value);              
			setResult("Spot_SC"+4+"_"+resultLabels[s]+"_Mean",0,value/count);
					
		}
		selectImage(matrix); close;
		updateResults();
		resultWindow=false;
	}

	selectWindow("Results"); 
	saveAs("Measurements",results);
}

function sortResults()
{
	resultLabels = getResultLabels();
	matrix = results2matrix(resultLabels);
	matrix2results(matrix,resultLabels,channels);
}

function getResultLabels()
{
	selectWindow("Results");
	ls 				= split(getInfo(),'\n');
	rr 				= split(ls[0],'\t'); 
	nparams 		= rr.length-1;			
	resultLabels 	= newArray(nparams);
	for(j=1;j<=nparams;j++){resultLabels[j-1]=rr[j];}
	return resultLabels;
}

function results2matrix(resultLabels)
{
	h = nResults;
	w = resultLabels.length;
	newImage("Matrix", "32-bit Black",w, h, 1);
	matrix = getImageID;
	for(j=0;j<w;j++)
	{
		for(r=0;r<h;r++)
		{
			v = getResult(resultLabels[j],r);
			selectImage(matrix);
			setPixel(j,r,v);
		}
	}
	run("Clear Results");
	return matrix;
}

function matrix2results(matrix,resultLabels,channels)
{
	selectImage(matrix);
	w = getWidth;
	h = getHeight;
	for(c=0;c<channels;c++)
	{
		start = c*h/channels;
		end = c*h/channels+h/channels;
		for(k=0;k<w;k++)
		{
			for(j=start;j<end;j++)
			{
				selectImage(matrix);
				p = getPixel(k,j);
				setResult(resultLabels[k]+"_MC"+c+1,j-start,p); // MC for measurement channel
			}
		}
	}
	selectImage(matrix); close;
	updateResults;
}

function toggleOverlay()
{	
	run("Select None"); 
	roiManager("deselect");
	roiManager("Show All without labels");
	if(Overlay.size == 0)run("From ROI Manager");
	else run("Remove Overlay");
}
	
function createOverlay(names)
{
	setForegroundColor(25, 25, 25);
	fields = names.length;
	print(fields,"images");
	for(i=0;i<fields;i++){
		prefix = names[i];
		file = prefix+suffix;
		setFileNames(prefix);
		print(i+1,"/",fields,":",prefix);
		path = dir+file;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		//open(path);
		id = getImageID;
		Stack.getDimensions(w,h,channels,slices,frames); 
		if(!Stack.isHyperStack && channels == 1){
			channels = slices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices=1 frames=1 display=Composite");
		}
		id = getImageID;
		if(spots_ch1){
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(spots_ch1_roi_set)){	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",spots_ch1_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
		
		if(spots_ch2){
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(spots_ch2_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",spots_ch2_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}

		if(spots_ch3){
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(spots_ch3_roi_set)){	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",spots_ch3_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}

		if(spots_ch4){
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(spots_ch4_roi_set)){	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",spots_ch4_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
	}
	run("Concatenate...", "all_open title=[Concatenated Stacks]");
	Stack.getDimensions(w,h,newchannels,slices,frames);
	for(c=1;c<=channels;c++){
		Stack.setChannel(c);
		Stack.setFrame(round(frames/2));
		resetMinAndMax;
	}
	range = pow(2,bitDepth);
	for(c=channels+1;c<=newchannels;c++){
		Stack.setChannel(c);
		setMinAndMax(0,range/2);
	}
	run("Make Composite");
}

function concatenateResults()
{
	prefixes  = scanFiles();
	index = 0;
	names = newArray(0);
	for(i = 0; i < prefixes.length; i++)
	{
		print(i+1,"/",prefixes.length,"...");
		results = output_dir+prefixes[i]+"_summary.txt";
		run("Results... ","open=["+results+"]");
		nr = nResults; print("...",nr,"results");
		for(r = index; r < index+nr; r++)names = Array.concat(names,prefixes[i]); 
			if(i > 0) {	
				labels = getResultLabels();
				matrix = results2matrix(labels);
				selectWindow("Results"); run("Close");
				IJ.renameResults("Summary","Results");

				//matrix2result
				selectImage(matrix);
				w = getWidth;
				h = getHeight;
				start = nResults;
				for(k = 0;k < w; k++)
				{
					for(j = 0; j < h; j++)
					{
						selectImage(matrix);
						p = getPixel(k,j);
						setResult(labels[k],j+start,p); 
					}
				}
				selectImage(matrix); close;
				updateResults;
			}
		selectWindow("Results"); 
		IJ.renameResults("Summary");
	}
	IJ.renameResults("Summary","Results");
	for(i = 0; i < names.length; i++)setResult("Label",i,names[i]);
	updateResults;
}