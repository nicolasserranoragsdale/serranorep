

import ij.IJ;
import ij.ImagePlus;
import ij.Undo;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.gui.OvalRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.filter.RankFilters;
import ij.process.Blitter;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;

public class My_Threshold_std_dev implements PlugIn {
	private double std_thresh;
	private int multfactor;
	private int bitDepth;
	private double k2;
	private double[] ContrastThresh;


	/** Ask for parameters and then execute.*/
	public void run(String arg) {
		// 1 - Obtain the currently active image:
		ImagePlus imp = IJ.getImage();

		if (null == imp){
			IJ.showMessage("There must be at least one image open");
			return;
		}
		bitDepth=imp.getBitDepth();
		std_thresh=5;
		multfactor=18;
		k2=6;

		int mean_value=128;
		int contrast_thresh=15;
		if (bitDepth == 16) {
			mean_value = 40000;
			std_thresh = 500;
			contrast_thresh=2500;
		}

		

		// 2 - Ask for parameters:
		GenericDialog gd = new GenericDialog("Auto Local Threshold with std_dev");
		String [] methods={"Bernsen","BernsenOtsu", "std_dev","stdwithBernsen",
				"BernsenKim","Otsu","BernsenContrast"};
		final Package p = getClass().getPackage();
		final String version = p == null ? null : p.getImplementationVersion();
		final String versionSuffix = version == null ? "" : " v" + version;
		gd.addMessage("Auto Local Threshold with std_dev" + versionSuffix);
		gd.addChoice("Method", methods, methods[0]);
		gd.addNumericField ("Radius",  5, 0);
		gd.addMessage ("Special parameters (if different from default)");
		gd.addNumericField ("Bernsen_contrast_threshold",  contrast_thresh, 0);
		gd.addNumericField ("Bernsen_mean_value",  mean_value, 0);
		if (imp.getStackSize()>1) {
			gd.addCheckbox("Stack",true);
		}

		gd.addNumericField("low std_dev threshold", std_thresh, 1);
		gd.addNumericField("multiplication factor for std_dev contrast threshold", multfactor, 0);
		gd.addNumericField("multiplication factor (Bersen Contrast)", k2, 0);

		gd.addMessage("Thresholded result is always shown in white [255].");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		// 3 - Retrieve parameters from the dialog
		String myMethod= gd.getNextChoice ();
		int radius = (int) gd.getNextNumber();
		double par1 = (double) gd.getNextNumber();
		double par2 = (double) gd.getNextNumber();
		boolean doIstack=false; 

		int stackSize = imp.getStackSize();
		if (stackSize>1)
			doIstack = gd.getNextBoolean ();

		std_thresh = gd.getNextNumber() ;
		multfactor=(int)gd.getNextNumber() ;
		k2=gd.getNextNumber() ;




		// 4 - Execute!
		long start;
		long end;
		
		if ((myMethod.equals("BernsenKim") || myMethod.equals("BernsenContrast")) && stackSize>150) {
			
			GenericDialog gd2 = new GenericDialog("Type of Contrast Threshold");
			gd2.addMessage("Same global threshold for all slices or individual threshold for each slice?");
			gd2.addCheckbox("global threshold?",true);
			gd2.showDialog();
			boolean global=gd2.getNextBoolean();
			
			start = System.currentTimeMillis()/1000;
			
			if (global) {
				
				if (myMethod.equals("BernsenKim")) {			
					ContrastThresh=new double[stackSize];
					
					if (bitDepth == 16) {
						for (int k=1; k<=stackSize; k++){
							imp.setSlice(k);
							ContrastThresh[k-1]=BernsenKim16(imp, radius, par1, par2);
						}
						double sum=0;
						for (int k=51; k<=stackSize-50; k++)
							sum+=ContrastThresh[k-1];
						double avgThresh=sum/((double)(stackSize-50-51+1));
						for (int k=1; k<=stackSize; k++) {
							imp.setSlice(k);
							Bernsen16(imp,  radius, avgThresh, par2);
						}
						imp.setSlice(1);
					}
						
					if (bitDepth == 8) {
						for (int k=1; k<=stackSize; k++){
							imp.setSlice(k);
							ContrastThresh[k-1]=BernsenKim(imp, radius, par1, par2);
						}
						double sum=0;
						for (int k=51; k<=stackSize-50; k++)
							sum+=ContrastThresh[k-1];
						double avgThresh=sum/((double)(stackSize-50-51+1));
						for (int k=1; k<=stackSize; k++) {
							imp.setSlice(k);
							Bernsen(imp,  radius, avgThresh, par2);
						}
						imp.setSlice(1);
					}
				}
				
				if (myMethod.equals("BernsenContrast")) {			
					ContrastThresh=new double[stackSize];
					
					if (bitDepth == 16) {
						for (int k=1; k<=stackSize; k++){
							imp.setSlice(k);
							ContrastThresh[k-1]=BernsenContrast16(imp, radius, par1, par2);
						}
						double sum=0;
						for (int k=51; k<=stackSize-50; k++)
							sum+=ContrastThresh[k-1];
						double avgThresh=sum/((double)(stackSize-50-51+1));
						for (int k=1; k<=stackSize; k++) {
							imp.setSlice(k);
							Bernsen16(imp,  radius, avgThresh, par2);
						}
						imp.setSlice(1);
					}
						
					if (bitDepth == 8) {
						for (int k=1; k<=stackSize; k++){
							imp.setSlice(k);
							ContrastThresh[k-1]=BernsenContrast(imp, radius, par1, par2);
						}
						double sum=0;
						for (int k=51; k<=stackSize-50; k++)
							sum+=ContrastThresh[k-1];
						double avgThresh=sum/((double)(stackSize-50-51+1));
						for (int k=1; k<=stackSize; k++) {
							imp.setSlice(k);
							Bernsen(imp,  radius, avgThresh, par2);
						}
						imp.setSlice(1);
					}
				}
			}
			
			else {
				
				if (myMethod.equals("BernsenKim")) {			
					double Contrast_Thresh;
					
					if (bitDepth == 16) {
						for (int k=1; k<=stackSize; k++){
							imp.setSlice(k);
							Contrast_Thresh=BernsenKim16(imp, radius, par1, par2);
							Bernsen16(imp,  radius, Contrast_Thresh, par2);
						}
						imp.setSlice(1);
					}
					if (bitDepth == 8) {
						for (int k=1; k<=stackSize; k++){
							imp.setSlice(k);
							Contrast_Thresh=BernsenKim(imp, radius, par1, par2);
							Bernsen(imp,  radius, Contrast_Thresh, par2);
						}
						imp.setSlice(1);
					}
				}
				
				if (myMethod.equals("BernsenContrast")) {			
					double Contrast_Thresh;
					
					if (bitDepth == 16) {
						for (int k=1; k<=stackSize; k++){
							imp.setSlice(k);
							Contrast_Thresh=BernsenContrast16(imp, radius, par1, par2);
							Bernsen16(imp,  radius, Contrast_Thresh, par2);
						}
						imp.setSlice(1);
					}
					if (bitDepth == 8) {
						for (int k=1; k<=stackSize; k++){
							imp.setSlice(k);
							Contrast_Thresh=BernsenContrast(imp, radius, par1, par2);
							Bernsen(imp,  radius, Contrast_Thresh, par2);
						}
						imp.setSlice(1);
					}
				}
			}
			
			end = System.currentTimeMillis()/1000;
		}


		
		

		else {
			start = System.currentTimeMillis()/1000;
			
			if (stackSize>1 && doIstack ) { 
				for (int k=1; k<=stackSize; k++){
					imp.setSlice(k);
					exec(imp, myMethod, radius, par1, par2);
				}
				
				imp.setSlice(1);
			}
			else { //just one slice
				exec(imp, myMethod, radius, par1, par2);
			}
			
			end = System.currentTimeMillis()/1000;
		}
		
		IJ.log("time(seconds): "+start+" "+end+"  "+(end-start));
		
	}
	


	
	public void exec(ImagePlus imp, String myMethod, int radius,  double par1, double par2) {

		// 0 - Check validity of parameters
		if (null == imp) return;
		ImageProcessor ip = imp.getProcessor();


		IJ.showStatus("Thresholding...");
		long startTime = System.currentTimeMillis();
		//1 Do it
		if (imp.getStackSize()==1){
			ip.snapshot();
			Undo.setup(Undo.FILTER, imp);
		}
		// Apply the selected algorithm
		if(myMethod.equals("Bernsen")){
			if (bitDepth == 16)
				Bernsen16(imp,  radius, par1, par2);
			if (bitDepth == 8)
				Bernsen(imp,  radius, par1, par2);
		}
		else if(myMethod.equals("BernsenOtsu")){
			BernsenOtsu(imp, radius, par1, par2);
		}
		else if(myMethod.equals("std_dev")){
			std(imp, radius);
		}
		else if(myMethod.equals("stdwithBernsen")){
			stdwithBernsen(imp, radius, par1, par2);
		}
		else if(myMethod.equals("BernsenKim")){
			if (bitDepth == 16) {
				double Contrast_Thresh=BernsenKim16(imp, radius, par1, par2);
				Bernsen16(imp,  radius, Contrast_Thresh, par2);
			
			}
			if (bitDepth == 8) {
				double Contrast_Thresh=BernsenKim(imp, radius, par1, par2);
				Bernsen(imp,  radius, Contrast_Thresh, par2);
			}
		}
		else if(myMethod.equals("BernsenContrast")){
			if (bitDepth == 16) {
				double Contrast_Thresh=BernsenContrast16(imp, radius, par1, par2);
				Bernsen16(imp,  radius, Contrast_Thresh, par2);
			
			}
			if (bitDepth == 8) {
				double Contrast_Thresh=BernsenContrast(imp, radius, par1, par2);
				Bernsen(imp,  radius, Contrast_Thresh, par2);
			}
		}
		else if(myMethod.equals("Otsu")){
			Otsu(imp, radius, par1, par2);
		}
		/*//IJ.showProgress((double)(255-i)/255);
		imp.updateAndDraw();
		imp.getProcessor().setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE);
		// 2 - Return the threshold and the image
		IJ.showStatus("\nDone " + (System.currentTimeMillis() - startTime) / 1000.0);
		return new Object[] {imp};*/
	}


	private void std(ImagePlus imp, int radius) {
		int w=imp.getWidth();
		int h=imp.getHeight();
		int position;
		int radiusx2=radius * 2;

		ImagePlus curr=imp;
		ImageConverter convrt=new ImageConverter(curr);
		convrt.convertToGray32();

		ImageProcessor ip=curr.getProcessor();
		float[] pixels = (float[]) ip.getPixels();
		int roiy;

		//
		float[] stds=new float[pixels.length];
		//

		Roi roi = new OvalRoi(0, 0, radiusx2, radiusx2);
		//ip.setRoi(roi);
		for (int y =0; y<h; y++){
			IJ.showProgress((double)(y)/(h-1)); // this method is slow, so let's show the progress bar
			roiy = y-radius;
			for (int x = 0; x<w; x++){
				roi.setLocation(x-radius,roiy);
				ip.setRoi(roi);
				//ip.setRoi(new OvalRoi(x-radius, roiy, radiusx2, radiusx2));
				position=x+y*w;


				float localstd=(float)ip.getStats().stdDev;
				stds[position]=localstd;

			}
		}
		for (position=0; position<w*h; position++) pixels[position]= stds[position]; //update with thresholded pixels
	}



	void stdwithBernsen(ImagePlus imp, int radius, double par1, double par2) {
		//only 8bit
		int w=imp.getWidth();
		int h=imp.getHeight();
		int position;
		int radiusx2=radius * 2;
		float object;
		float backg;
		float mean_contrast=128;

		if (par2!=128) {
			IJ.log("Bernsen: changed mean_contrast from :"+ mean_contrast + "  to:" + par2);
			mean_contrast= (float) par2;
		}

		object =  (float) 0xff;
		backg = (float) 0;

		ImagePlus curr=imp;
		ImageConverter convrt=new ImageConverter(curr);
		convrt.convertToGray32();

		ImageProcessor ip=curr.getProcessor();
		float[] pixels = (float[]) ip.getPixels();
		int roiy;
		float[] copy=new float[pixels.length];

		Roi roi = new OvalRoi(0, 0, radiusx2, radiusx2);
		//ip.setRoi(roi);
		for (int y =0; y<h; y++){
			IJ.showProgress((double)(y)/(h-1)); // this method is slow, so let's show the progress bar
			roiy = y-radius;
			for (int x = 0; x<w; x++){
				roi.setLocation(x-radius,roiy);
				ip.setRoi(roi);
				//ip.setRoi(new OvalRoi(x-radius, roiy, radiusx2, radiusx2));
				position=x+y*w;


				float localstd=(float)ip.getStats().stdDev;

				float localmax=(float) ip.getStats().max;
				float localmin=(float) ip.getStats().min;
				float midgray=(localmax+localmin)/2;

				//standard deviation will become local threshold
				if ( localstd < std_thresh )
					copy[position] = ( midgray >= mean_contrast ) ? object :  backg;  //Low contrast region
				else
					copy[position] = (pixels[position] > midgray ) ? object : backg;
			}
		}
		for (position=0; position<w*h; position++) pixels[position]=copy[position];
	}

	//only 8bit
	void BernsenOtsu(ImagePlus imp, int radius, double par1, double par2) {
		int[] data;
		int w=imp.getWidth();
		int h=imp.getHeight();
		int position;
		int radiusx2=radius * 2;
		ImageProcessor ip=imp.getProcessor();
		byte[] pixels = (byte []) ip.getPixels();
		byte[] pixelsOut = new byte[pixels.length]; // need this to avoid changing the image data (and further histograms)
		byte object;
		byte backg;
		float contrast_threshold=15;
		float mean_contrast=128;

		if (par1!=15) {
			IJ.log("Bernsen: changed contrast_threshold from :"+ contrast_threshold + "  to:" + par1);
			contrast_threshold= (float)par1;
		}

		if (par2!=128) {
			IJ.log("Bernsen: changed mean_contrast from :"+ mean_contrast + "  to:" + par2);
			mean_contrast= (float) par2;
		}


		object =  (byte) 0xff;
		backg =   (byte) 0;



		int k,kStar;  // k = the current threshold; kStar = optimal threshold
		int N1, N;    // N1 = # points with intensity <=k; N = total number of points
		double BCV, BCVmax; // The current Between Class Variance and maximum BCV
		double num, denom;  // temporary bookeeping
		int Sk;  // The total intensity for all histogram points <=k
		int S, L=256; // The total intensity of the image. Need to hange here if modifying for >8 bits images
		int roiy;

		Roi roi = new OvalRoi(0, 0, radiusx2, radiusx2);
		//ip.setRoi(roi);
		for (int y =0; y<h; y++){
			IJ.showProgress((double)(y)/(h-1)); // this method is slow, so let's show the progress bar
			roiy = y-radius;
			for (int x = 0; x<w; x++){
				roi.setLocation(x-radius,roiy);
				ip.setRoi(roi);
				//ip.setRoi(new OvalRoi(x-radius, roiy, radiusx2, radiusx2));
				position=x+y*w;
				data = ip.getHistogram();

				double localmax=ip.getStats().max;
				double localmin=ip.getStats().min;
				double local_contrast=localmax-localmin;
				double midgray=(localmax+localmin)/2;

				// Initialize values:
				S = N = 0;
				for (k=0; k<L; k++) {
					S += k * data[k];	// Total histogram intensity
					N += data[k];		// Total number of data points
				}

				Sk = 0;
				N1 = data[0]; // The entry for zero intensity
				BCV = 0;
				BCVmax=0;
				kStar = 0;

				// Look at each possible threshold value,
				// calculate the between-class variance, and decide if it's a max
				for (k=1; k<L-1; k++) { // No need to check endpoints k = 0 or k = L-1
					Sk += k * data[k];
					N1 += data[k];

					// The float casting here is to avoid compiler warning about loss of precision and
					// will prevent overflow in the case of large saturated images
					denom = (double)( N1) * (N - N1); // Maximum value of denom is (N^2)/4 =  approx. 3E10

					if (denom != 0 ){
						// Float here is to avoid loss of precision when dividing
						num = ( (double)N1 / N ) * S - Sk; 	// Maximum value of num =  255*N = approx 8E7
						BCV = (num * num) / denom;
					}
					else
						BCV = 0;

					if (BCV >= BCVmax){ // Assign the best threshold found so far
						BCVmax = BCV;
						kStar = k;
					}
				}
				// kStar += 1;	// Use QTI convention that intensity -> 1 if intensity >= k
				// (the algorithm was developed for I-> 1 if I <= k.)
				//return kStar;

				if (local_contrast < contrast_threshold)
					pixelsOut[position] = ( midgray >= mean_contrast ) ? object :  backg;
				else	
					pixelsOut[position] = ((int) (pixels[position]&0xff)>kStar) ? object : backg;
			}
		}
		for (position=0; position<w*h; position++) pixels[position]=pixelsOut[position]; //update with thresholded pixels

	}



	double BernsenKim(ImagePlus imp, int radius, double par1, double par2) {
		ImageProcessor ip=imp.getProcessor();
		ImageProcessor ipvar=ip.duplicate().convertToFloat();
		//ImagePlus curr=imp.duplicate();
		//ImageConverter convrt=new ImageConverter(curr);
		//convrt.convertToGray32();

		//ImagePlus varimp=curr;
		//ImageProcessor ipvar=varimp.getProcessor();
		RankFilters rf=new RankFilters();
		rf.rank(ipvar, radius, rf.VARIANCE); 

		byte[] pixels = (byte [])ip.getPixels();
		float[] var = (float [])ipvar.getPixels();


		//
		ImageProcessor ipMax=ip.duplicate();
		rf.rank(ipMax, radius, rf.MAX);// Maximum

		ImageProcessor ipMin=ip.duplicate();
		rf.rank(ipMin, radius, rf.MIN); //Minimum


		byte[] max = (byte [])ipMax.getPixels();
		byte[] min = (byte [])ipMin.getPixels();
		//

		int size=pixels.length;
		int numlowstds=0;
		double sumlowstds=0;

		for(int i=0;i<size;i++) {
			double localstd=Math.sqrt(var[i]);
			int local_contrast = (int)((max[i]&0xff) -(min[i]&0xff));
			if (local_contrast < par1 && localstd!=0) {
				numlowstds++;
				sumlowstds=sumlowstds+localstd;
			}
		}

		double meanlowstds=sumlowstds/(double)numlowstds;
		double local_threshold=multfactor*meanlowstds;
		return local_threshold;
		//Bernsen(imp,  radius, local_threshold, par2);
	}




	double BernsenKim16(ImagePlus imp, int radius, double par1, double par2) {
		ImageProcessor ip=imp.getProcessor();
		ImageProcessor ipvar=ip.duplicate().convertToFloat();
		//ImagePlus curr=imp.duplicate();
		//ImageConverter convrt=new ImageConverter(curr);
		//convrt.convertToGray32();

		//ImagePlus varimp=curr;
		//ImageProcessor ipvar=varimp.getProcessor();
		RankFilters rf=new RankFilters();
		rf.rank(ipvar, radius, rf.VARIANCE); 

		short[] pixels = (short [])ip.getPixels();
		float[] var = (float [])ipvar.getPixels();

		ImageProcessor ipMax=ip.duplicate();
		rf.rank(ipMax, radius, rf.MAX);// Maximum

		ImageProcessor ipMin=ip.duplicate();
		rf.rank(ipMin, radius, rf.MIN); //Minimum


		short[] max = (short [])ipMax.getPixels();
		short[] min = (short [])ipMin.getPixels();

		int size=pixels.length;
		int numlowstds=0;
		double sumlowstds=0;

		for(int i=0;i<size;i++) {
			double localstd=Math.sqrt(var[i]);
			int local_contrast = (int)((max[i]&0xffff) -(min[i]&0xffff));
			if (local_contrast < par1 && localstd!=0) {
				numlowstds++;
				sumlowstds=sumlowstds+localstd;
			}
		}

		double meanlowstds=sumlowstds/(double)numlowstds;
		double local_threshold=multfactor*meanlowstds;
		
		return local_threshold;
		//Bernsen16(imp,  radius, local_threshold, par2);
	}

	void Otsu(ImagePlus imp, int radius,  double par1, double par2) {
		// Otsu's threshold algorithm
		// C++ code by Jordan Bevik <Jordan.Bevic@qtiworld.com>
		// ported to ImageJ plugin by G.Landini. Same algorithm as in Auto_Threshold, this time on local circular regions
		int[] data;
		int w=imp.getWidth();
		int h=imp.getHeight();
		int position;
		int radiusx2=radius * 2;
		ImageProcessor ip=imp.getProcessor();
		byte[] pixels = (byte []) ip.getPixels();
		byte[] pixelsOut = new byte[pixels.length]; // need this to avoid changing the image data (and further histograms)
		byte object;
		byte backg;


		object =  (byte) 0xff;
		backg =   (byte) 0;


		int k,kStar;  // k = the current threshold; kStar = optimal threshold
		int N1, N;    // N1 = # points with intensity <=k; N = total number of points
		double BCV, BCVmax; // The current Between Class Variance and maximum BCV
		double num, denom;  // temporary bookeeping
		int Sk;  // The total intensity for all histogram points <=k
		int S, L=256; // The total intensity of the image. Need to hange here if modifying for >8 bits images
		int roiy;

		Roi roi = new OvalRoi(0, 0, radiusx2, radiusx2);
		//ip.setRoi(roi);
		for (int y =0; y<h; y++){
			IJ.showProgress((double)(y)/(h-1)); // this method is slow, so let's show the progress bar
			roiy = y-radius;
			for (int x = 0; x<w; x++){
				roi.setLocation(x-radius,roiy);
				ip.setRoi(roi);
				//ip.setRoi(new OvalRoi(x-radius, roiy, radiusx2, radiusx2));
				position=x+y*w;
				data = ip.getHistogram();

				// Initialize values:
				S = N = 0;
				for (k=0; k<L; k++){
					S += k * data[k];	// Total histogram intensity
					N += data[k];		// Total number of data points
				}

				Sk = 0;
				N1 = data[0]; // The entry for zero intensity
				BCV = 0;
				BCVmax=0;
				kStar = 0;

				// Look at each possible threshold value,
				// calculate the between-class variance, and decide if it's a max
				for (k=1; k<L-1; k++) { // No need to check endpoints k = 0 or k = L-1
					Sk += k * data[k];
					N1 += data[k];

					// The float casting here is to avoid compiler warning about loss of precision and
					// will prevent overflow in the case of large saturated images
					denom = (double)( N1) * (N - N1); // Maximum value of denom is (N^2)/4 =  approx. 3E10

					if (denom != 0 ){
						// Float here is to avoid loss of precision when dividing
						num = ( (double)N1 / N ) * S - Sk; 	// Maximum value of num =  255*N = approx 8E7
						BCV = (num * num) / denom;
					}
					else
						BCV = 0;

					if (BCV >= BCVmax){ // Assign the best threshold found so far
						BCVmax = BCV;
						kStar = k;
					}
				}
				// kStar += 1;	// Use QTI convention that intensity -> 1 if intensity >= k
				// (the algorithm was developed for I-> 1 if I <= k.)
				//return kStar;
				pixelsOut[position] = ((int) (pixels[position]&0xff)>kStar) ? object : backg;
			}
		}
		for (position=0; position<w*h; position++) pixels[position]=pixelsOut[position]; //update with thresholded pixels
	}

	void Bernsen(ImagePlus imp, int radius,  double par1, double par2) {
		// Bernsen recommends WIN_SIZE = 31 and CONTRAST_THRESHOLD = 15.
		//  1) Bernsen J. (1986) "Dynamic Thresholding of Grey-Level Images" 
		//    Proc. of the 8th Int. Conf. on Pattern Recognition, pp. 1251-1255
		//  2) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
		//   Techniques and Quantitative Performance Evaluation" Journal of 
		//   Electronic Imaging, 13(1): 146-165 
		//  http://citeseer.ist.psu.edu/sezgin04survey.html
		// Ported to ImageJ plugin from E Celebi's fourier_0.8 routines
		// This version uses a circular local window, instead of a rectagular one
		ImagePlus Maximp, Minimp;
		ImageProcessor ip=imp.getProcessor(), ipMax, ipMin;
		float contrast_threshold=15;
		int local_contrast;
		int mid_gray;
		byte object;
		byte backg;
		int temp;
		float mean_value=128;

		if (par1!=15) {
			//IJ.log("Bernsen: changed contrast_threshold from :"+ contrast_threshold + "  to:" + par1);
			IJ.log(""+ par1);
			contrast_threshold= (float) par1;
		}

		if (par2!=128) {
			//IJ.log("Bernsen: changed mean_value from :"+ mean_value + "  to:" + par2);
			mean_value= (float) par2;
		}


		object =  (byte) 0xff;
		backg =   (byte) 0;



		Maximp=duplicateImage(ip);
		ipMax=Maximp.getProcessor();
		RankFilters rf=new RankFilters();
		rf.rank(ipMax, radius, rf.MAX);// Maximum
		//Maximp.show();
		Minimp=duplicateImage(ip);
		ipMin=Minimp.getProcessor();
		rf.rank(ipMin, radius, rf.MIN); //Minimum
		//Minimp.show();
		byte[] pixels = (byte [])ip.getPixels();
		byte[] max = (byte [])ipMax.getPixels();
		byte[] min = (byte [])ipMin.getPixels();

		int porecount=0;
		int whitecount=0;	
		for (int i=0; i<pixels.length; i++) {
			local_contrast = (int)((max[i]&0xff) -(min[i]&0xff));
			mid_gray =(int) ((min[i]&0xff) + (max[i]&0xff) )/ 2;
			temp=(int) (pixels[i] & 0x0000ff);

			if ( local_contrast < contrast_threshold ) { //Low contrast region
				if (mid_gray >= mean_value && temp!=0) {
					whitecount++;
					pixels[i]=object;
				}
				if (mid_gray < mean_value && temp!=0) {
					porecount++;
					pixels[i]=backg;
				}
				//pixels[i] = ( mid_gray >= mean_value ) ? object :  backg;
			}
			else {
				if (temp >= mid_gray && temp!=0) {
					whitecount++;
					pixels[i]=object;
				}
				if (temp < mid_gray && temp!=0) {
					porecount++;
					pixels[i]=backg;
				}
				//pixels[i] = (temp >= mid_gray ) ? object : backg;
			}
		}   
		//IJ.log(""+(float)(porecount)*100/(float)(whitecount+porecount));
		//imp.updateAndDraw();
		return;
	}


	void Bernsen16(ImagePlus imp, int radius,  double par1, double par2) {
		ImageProcessor ip=imp.getProcessor(), ipMax, ipMin;
		float contrast_threshold=2500;
		int local_contrast;
		float mid_gray;
		byte object;
		byte backg;
		int temp;
		float mean_value=40000;

		if (par1!=2500) {
			//IJ.log("Bernsen: changed contrast_threshold from :"+ contrast_threshold + "  to:" + par1);
			IJ.log(""+ par1);
			contrast_threshold= (float) par1;
		}

		if (par2!=40000) {
			//IJ.log("Bernsen: changed mean_value from :"+ mean_value + "  to:" + par2);
			mean_value= (float) par2;
		}


		object =  (byte) 0xff;
		backg =   (byte) 0;



		//Maximp=imp.duplicate();
		//ipMax=Maximp.getProcessor();
		ipMax=ip.duplicate();
		RankFilters rf=new RankFilters();
		rf.rank(ipMax, radius, rf.MAX);// Maximum
		//Maximp.show();

		//Minimp=imp.duplicate();
		//ipMin=Minimp.getProcessor();
		ipMin=ip.duplicate();
		rf.rank(ipMin, radius, rf.MIN); //Minimum
		//Minimp.show();
		short[] pixels = (short [])ip.getPixels();
		short[] max = (short [])ipMax.getPixels();
		short[] min = (short [])ipMin.getPixels();

		int porecount=0;
		int whitecount=0;

		for (int i=0; i<pixels.length; i++) {
			local_contrast = ((max[i]& 0xffff) -(min[i])& 0xffff);
			mid_gray =(float) (((min[i]& 0xffff) + (max[i]& 0xffff) )/ 2.0);
			temp=(int) (pixels[i] & 0xffff);

			if ( local_contrast < contrast_threshold ) { //Low contrast region
				if (mid_gray >= mean_value && temp!=0) {
					whitecount++;
					pixels[i]=object;
				}
				if (mid_gray < mean_value && temp!=0) {
					porecount++;
					pixels[i]=backg;
				}
				//pixels[i] = ( mid_gray >= mean_value ) ? object :  backg;
			}
			else {
				if (temp >= mid_gray && temp!=0) {
					whitecount++;
					pixels[i]=object;
				}
				if (temp < mid_gray && temp!=0) {
					porecount++;
					pixels[i]=backg;
				}
				//pixels[i] = (temp >= mid_gray ) ? object : backg;
			}
		}
		//IJ.log(""+(float)(porecount)*100/(float)(whitecount+porecount));
		//imp.updateAndDraw();
		return;
	}

	double BernsenContrast16(ImagePlus imp, int radius, double par1, double par2) {
		ImageProcessor ip=imp.getProcessor(), ipMax, ipMin;
		int local_contrast;

		ipMax=ip.duplicate();
		RankFilters rf=new RankFilters();
		rf.rank(ipMax, radius, rf.MAX);

		ipMin=ip.duplicate();
		rf.rank(ipMin, radius, rf.MIN);

		short[] pixels = (short [])ip.getPixels();
		short[] max = (short [])ipMax.getPixels();
		short[] min = (short [])ipMin.getPixels();

		double sum_lowcon=0;
		int num_lowcon=0;

		for (int i=0; i<pixels.length; i++) {
			local_contrast = ((max[i]& 0xffff) -(min[i]& 0xffff));

			if ( 1 <local_contrast && local_contrast< par1) {
				sum_lowcon=sum_lowcon+local_contrast;
				num_lowcon++;
			}
		}  
		double mean_lowcon=sum_lowcon/(double)num_lowcon;
		//IJ.log("mean contrast of low contrast regions:"+mean_lowcon);

		double sumsquares_lowcon=0;
		for (int i=0; i<pixels.length; i++) {
			local_contrast = ((max[i]& 0xffff) -(min[i]& 0xffff));

			if ( 1 <local_contrast && local_contrast< par1)
				sumsquares_lowcon+=(local_contrast-mean_lowcon)*(local_contrast-mean_lowcon);
		}
		double std_lowcon=Math.sqrt(sumsquares_lowcon/((double)num_lowcon-1));

		//IJ.log("std:"+std_lowcon);

		double contrast_thresh=mean_lowcon+k2*std_lowcon;
		
		return contrast_thresh;
		//Bernsen16(imp,  radius, contrast_thresh, par2);
	}


	double BernsenContrast(ImagePlus imp, int radius, double par1, double par2) {
		ImageProcessor ip=imp.getProcessor(), ipMax, ipMin;
		int local_contrast;

		ipMax=ip.duplicate();
		RankFilters rf=new RankFilters();
		rf.rank(ipMax, radius, rf.MAX);

		ipMin=ip.duplicate();
		rf.rank(ipMin, radius, rf.MIN);

		byte[] pixels = (byte [])ip.getPixels();
		byte[] max = (byte [])ipMax.getPixels();
		byte[] min = (byte [])ipMin.getPixels();

		double sum_lowcon=0;
		int num_lowcon=0;

		for (int i=0; i<pixels.length; i++) {
			local_contrast = ((max[i]& 0xff) -(min[i]& 0xff));

			if ( 1 <local_contrast && local_contrast< par1) {
				sum_lowcon=sum_lowcon+local_contrast;
				num_lowcon++;
			}
		}  
		double mean_lowcon=sum_lowcon/(double)num_lowcon;
		//IJ.log("mean contrast of low contrast regions:"+mean_lowcon);

		double sumsquares_lowcon=0;
		for (int i=0; i<pixels.length; i++) {
			local_contrast = ((max[i]& 0xff) -(min[i]& 0xff));

			if ( 1 <local_contrast && local_contrast< par1)
				sumsquares_lowcon+=(local_contrast-mean_lowcon)*(local_contrast-mean_lowcon);
		}
		double std_lowcon=Math.sqrt(sumsquares_lowcon/((double)num_lowcon-1));

		//IJ.log("std:"+std_lowcon);

		double contrast_thresh=mean_lowcon+k2*std_lowcon;
		
		return contrast_thresh;
		//Bernsen(imp,  radius, contrast_thresh, par2);
	}


	//experimental. Need to remove outliers.
	/*void FindThreshold(ImagePlus imp, int radius, double par1, double par2) {
		ImageProcessor ip=imp.getProcessor();
		short[] pixels = (short [])ip.getPixels();
		int[] pixelsconvrt=new int[pixels.length];
		for (int i=0;i<pixels.length;i++)
			pixelsconvrt[i]=(pixels[i] & 0xffff);

		Arrays.sort(pixelsconvrt);

		int maxdiff=0;
		double bestmidgray=0;
		int value=0;
		for (int i=0;i<pixelsconvrt.length-1;i++) {
			int localdiff=pixelsconvrt[i+1]-pixelsconvrt[i];
			boolean nonzero=true;
			if (pixelsconvrt[i+1] < 1 || pixelsconvrt[i] < 1)
				nonzero=false;
			if (localdiff > maxdiff && nonzero) {
				maxdiff=localdiff;
				bestmidgray=(double)(pixelsconvrt[i+1]+pixelsconvrt[i])/2;
				value=pixelsconvrt[i];
			}
		}

		IJ.log("maximum difference: "+maxdiff);
		IJ.log("midgray: "+bestmidgray);
		IJ.log("value: "+value);
	}*/

	private ImagePlus duplicateImage(ImageProcessor iProcessor){
		int w=iProcessor.getWidth();
		int h=iProcessor.getHeight();
		ImagePlus iPlus=NewImage.createByteImage("Image", w, h, 1, NewImage.FILL_BLACK);
		ImageProcessor imageProcessor=iPlus.getProcessor();
		imageProcessor.copyBits(iProcessor, 0,0, Blitter.COPY);
		return iPlus;
	} 

}
