import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Undo;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.gui.OvalRoi;
import ij.gui.Roi;
import ij.gui.YesNoCancelDialog;
import ij.plugin.CanvasResizer;
import ij.plugin.ContrastEnhancer;
import ij.plugin.MontageMaker;
import ij.plugin.PlugIn;
import ij.plugin.filter.RankFilters;
import ij.process.Blitter;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;

public class My_Threshold_std_dev implements PlugIn {
	private double std_thresh;
	private int multfactor;
	private int bitDepth;
	
	
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
		
		int mean_contrast=128;
		int contrast_thresh=15;
		if (bitDepth == 16) {
			mean_contrast = 40000;
			std_thresh = 500;
			contrast_thresh=2500;
		}

		//if (imp.getBitDepth()!=8) {
		//	IJ.showMessage("Error", "Only 8-bit images are supported");
		//	return;
		//}

		// 2 - Ask for parameters:
		GenericDialog gd = new GenericDialog("Auto Local Threshold with std_dev");
		String [] methods={"Bernsen","BernsenOtsu", "std_dev","stdwithBernsen",
				"stdwithBernsenKim","Otsu"};
		final Package p = getClass().getPackage();
		final String version = p == null ? null : p.getImplementationVersion();
		final String versionSuffix = version == null ? "" : " v" + version;
		gd.addMessage("Auto Local Threshold with std_dev" + versionSuffix);
		gd.addChoice("Method", methods, methods[0]);
		gd.addNumericField ("Radius",  5, 0);
		gd.addMessage ("Special parameters (if different from default)");
		gd.addNumericField ("Bernsen_contrast_threshold",  contrast_thresh, 0);
		gd.addNumericField ("Bernsen_mean_contrast",  mean_contrast, 0);
		if (imp.getStackSize()>1) {
			gd.addCheckbox("Stack",true);
		}
		
		gd.addNumericField("low std_dev threshold", std_thresh, 1);
		gd.addNumericField("multiplication factor for std_dev contrast threshold", multfactor, 0);
		
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
		
		// 4 - Execute!
		//long start = System.currentTimeMillis();

		if (stackSize>1 &&  doIstack ) { //whole stack
			//				if (doIstackHistogram) {// one global histogram
			//					Object[] result = exec(imp, myMethod, noWhite, noBlack, doIwhite, doIset, doIlog, doIstackHistogram );
			//				}
			//				else{ // slice by slice
			for (int k=1; k<=stackSize; k++){
				imp.setSlice(k);
				Object[] result = exec(imp, myMethod, radius, par1, par2);
			}
			//				}
			imp.setSlice(1);
		}
		else { //just one slice
			Object[] result = exec(imp, myMethod, radius, par1, par2);
		}
		// 5 - If all went well, show the image:
		// not needed here as the source image is binarised 
	}
	//IJ.showStatus(IJ.d2s((System.currentTimeMillis()-start)/1000.0, 2)+" seconds");


	/** Execute the plugin functionality: duplicate and scale the given image.
	 * @return an Object[] array with the name and the scaled ImagePlus.
	 * Does NOT show the new, image; just returns it. */
	public Object[] exec(ImagePlus imp, String myMethod, int radius,  double par1, double par2) {

		// 0 - Check validity of parameters
		if (null == imp) return null;
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
			std(imp, radius, par1, par2);
		}
		else if(myMethod.equals("stdwithBernsen")){
			stdwithBernsen(imp, radius, par1, par2);
		}
		else if(myMethod.equals("stdwithBernsenKim")){
			if (bitDepth == 16)
				BernsenKim16(imp, radius, par1, par2);
			if (bitDepth == 8)
				BernsenKim(imp, radius, par1, par2);
		}
		else if(myMethod.equals("Otsu")){
			Otsu(imp, radius, par1, par2);
		}
		//IJ.showProgress((double)(255-i)/255);
		imp.updateAndDraw();
		imp.getProcessor().setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE);
		// 2 - Return the threshold and the image
		IJ.showStatus("\nDone " + (System.currentTimeMillis() - startTime) / 1000.0);
		return new Object[] {imp};
	}


	private void std(ImagePlus imp, int radius, double par1, double par2) {
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
	
	

	void BernsenKim(ImagePlus imp, int radius, double par1, double par2) {
		int w=imp.getWidth();
		int h=imp.getHeight();
		int radiusx2=radius * 2;

		
		//ImagePlus curr=imp;
		//ImageConverter convrt=new ImageConverter(curr);
		//convrt.convertToGray32();

		ImageProcessor ip=imp.getProcessor();
		int roiy;
		double sumlowstds=0;
		int numlowstds=0;

		Roi roi = new OvalRoi(0, 0, radiusx2, radiusx2);
		//ip.setRoi(roi);
		for (int y =0; y<h; y++){
			IJ.showProgress((double)(y)/(h-1)); // this method is slow, so let's show the progress bar
			roiy = y-radius;
			for (int x = 0; x<w; x++){
				roi.setLocation(x-radius,roiy);
				ip.setRoi(roi);
				//ip.setRoi(new OvalRoi(x-radius, roiy, radiusx2, radiusx2));

				double localstd=(double)ip.getStats().stdDev;
				if (localstd < std_thresh) {
					numlowstds++;
					sumlowstds=sumlowstds+localstd;
				}
				
			}
		}
		double meanlowstds=sumlowstds/(double)numlowstds;
		double local_threshold=multfactor*meanlowstds;
		
		Bernsen(imp,  radius, local_threshold, par2);
	}
	
	
	
	
	void BernsenKim16(ImagePlus imp, int radius, double par1, double par2) {
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
		
		int size=pixels.length;
		int numlowstds=0;
		double sumlowstds=0;
		
		for(int i=0;i<size;i++) {
			double localstd=Math.sqrt(var[i]);
			if (localstd < std_thresh) {
				numlowstds++;
				sumlowstds=sumlowstds+localstd;
			}
		}
		
		double meanlowstds=sumlowstds/(double)numlowstds;
		double local_threshold=multfactor*meanlowstds;
		
		Bernsen16(imp,  radius, local_threshold, par2);
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
		float mean_contrast=128;
		
		if (par1!=15) {
			IJ.log("Bernsen: changed contrast_threshold from :"+ contrast_threshold + "  to:" + par1);
			contrast_threshold= (float) par1;
		}
		
		if (par2!=128) {
			IJ.log("Bernsen: changed mean_contrast from :"+ mean_contrast + "  to:" + par2);
			mean_contrast= (float) par2;
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

		for (int i=0; i<pixels.length; i++) {
			local_contrast = (int)((max[i]&0xff) -(min[i]&0xff));
			mid_gray =(int) ((min[i]&0xff) + (max[i]&0xff) )/ 2;
			temp=(int) (pixels[i] & 0x0000ff);
			if ( local_contrast < contrast_threshold )
				pixels[i] = ( mid_gray >= mean_contrast ) ? object :  backg;  //Low contrast region
			else
				pixels[i] = (temp >= mid_gray ) ? object : backg;
		}    
		//imp.updateAndDraw();
		return;
	}
	
	
	void Bernsen16(ImagePlus imp, int radius,  double par1, double par2) {
		ImagePlus Maximp, Minimp;
		ImageProcessor ip=imp.getProcessor(), ipMax, ipMin;
		float contrast_threshold=2500;
		int local_contrast;
		float mid_gray;
		byte object;
		byte backg;
		int temp;
		float mean_contrast=40000;
		
		if (par1!=2500) {
			IJ.log("Bernsen: changed contrast_threshold from :"+ contrast_threshold + "  to:" + par1);
			contrast_threshold= (float) par1;
		}
		
		if (par2!=40000) {
			IJ.log("Bernsen: changed mean_contrast from :"+ mean_contrast + "  to:" + par2);
			mean_contrast= (float) par2;
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

		for (int i=0; i<pixels.length; i++) {
			local_contrast = ((max[i]& 0xffff) -(min[i])& 0xffff);
			mid_gray =(float) (((min[i]& 0xffff) + (max[i]& 0xffff) )/ 2.0);
			temp=(int) (pixels[i] & 0xffff);
			if ( local_contrast < contrast_threshold )
				pixels[i] = ( mid_gray >= mean_contrast ) ? object :  backg;  //Low contrast region
			else
				pixels[i] = (temp >= mid_gray ) ? object : backg;
		}    
		//imp.updateAndDraw();
		return;
	}
	
	
	void Bernsen16og(ImagePlus imp, int radius,  double par1, double par2) {
		float contrast_threshold=2500;
		byte object;
		byte backg;
		float mean_contrast=40000;
		
		int position;
		int radiusx2=radius * 2;
		int w=imp.getWidth();
		int h=imp.getHeight();
		
		if (par1!=2500) {
			IJ.log("Bernsen: changed contrast_threshold from :"+ contrast_threshold + "  to:" + par1);
			contrast_threshold= (float) par1;
		}
		
		if (par2!=40000) {
			IJ.log("Bernsen: changed mean_contrast from :"+ mean_contrast + "  to:" + par2);
			mean_contrast= (float) par2;
		}
		

		object =  (byte) 0xff;
		backg =   (byte) 0;

		ImageProcessor ip=imp.getProcessor();
		short[] pixels =  (short[]) ip.getPixels();
		int roiy;
		byte[] copy=new byte[pixels.length];

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


				float localmax=(float) ip.getStats().max;
				float localmin=(float) ip.getStats().min;
				float local_contrast=localmax-localmin;
				float midgray=(localmax+localmin)/2;
				
				if ( local_contrast < contrast_threshold )
					copy[position] = ( midgray > mean_contrast ) ? object :  backg;  //Low contrast region
				else
					copy[position] = ((pixels[position]& 0xffff) >= midgray ) ? object : backg;
			}
		}
		for (position=0; position<w*h; position++) pixels[position]= copy[position];
	}
	


	private ImagePlus duplicateImage(ImageProcessor iProcessor){
		int w=iProcessor.getWidth();
		int h=iProcessor.getHeight();
		ImagePlus iPlus=NewImage.createByteImage("Image", w, h, 1, NewImage.FILL_BLACK);
		ImageProcessor imageProcessor=iPlus.getProcessor();
		imageProcessor.copyBits(iProcessor, 0,0, Blitter.COPY);
		return iPlus;
	} 

}
