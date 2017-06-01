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

// AutoLocalThreshold segmentation 
// Following the guidelines at http://pacific.mpi-cbg.de/wiki/index.php/PlugIn_Design_Guidelines
// ImageJ plugin by G. Landini at bham. ac. uk
// 1.0  15/Apr/2009
// 1.1  01/Jun/2009
// 1.2  25/May/2010
// 1.3  1/Nov/2011 added constant offset to Niblack's method (request)
// 1.4  2/Nov/2011 Niblack's new constant should be subtracted to match mean,mode and midgrey methods. Midgrey method had the wrong constant sign.
// 1.5  18/Nov/2013 added 3 new local thresholding methdos: Constrast, Otsu and Phansalkar
// 1.6  18/Feb/2016 Stefan Helfrich fixed a typo. 
// 1.7  21/Jun/2016 Arttu Miettinen found that the the standard deviation in the Phansalkar method was not being computed properly

public class std_dev implements PlugIn {
	/** Ask for parameters and then execute.*/
	public void run(String arg) {
		// 1 - Obtain the currently active image:
		ImagePlus imp = IJ.getImage();

		if (null == imp){
			IJ.showMessage("There must be at least one image open");
			return;
		}

		if (imp.getBitDepth()!=8) {
			IJ.showMessage("Error", "Only 8-bit images are supported");
			return;
		}

		// 2 - Ask for parameters:
		GenericDialog gd = new GenericDialog("Auto Local Threshold");
		String [] methods={"Try all", "Bernsen","BernsenOtsu", "std_dev","stdwithBernsen",
				"stdwithBernsenKim"};
		final Package p = getClass().getPackage();
		final String version = p == null ? null : p.getImplementationVersion();
		final String versionSuffix = version == null ? "" : " v" + version;
		gd.addMessage("Auto Local Threshold" + versionSuffix);
		gd.addChoice("Method", methods, methods[0]);
		gd.addNumericField ("Radius",  15, 0);
		gd.addMessage ("Special parameters (if different from default)");
		gd.addNumericField ("Parameter_1",  0, 0);
		gd.addNumericField ("Parameter_2",  0, 0);
		gd.addCheckbox("White objects on black background",true);
		if (imp.getStackSize()>1) {
			gd.addCheckbox("Stack",false);
		}
		gd.addMessage("Thresholded result is always shown in white [255].");
		gd.showDialog();
		if (gd.wasCanceled()) return;

		// 3 - Retrieve parameters from the dialog
		String myMethod= gd.getNextChoice ();
		int radius = (int) gd.getNextNumber();
		double par1 = (double) gd.getNextNumber();
		double par2 = (double) gd.getNextNumber();
		boolean doIwhite = gd.getNextBoolean ();
		boolean doIstack=false; 

		int stackSize = imp.getStackSize();
		if (stackSize>1)
			doIstack = gd.getNextBoolean ();

		// 4 - Execute!
		//long start = System.currentTimeMillis();
		if(myMethod.equals("Try all")){
			ImageProcessor ip = imp.getProcessor();
			int xe = ip.getWidth();
			int ye = ip.getHeight();
			int ml = methods.length;
			ImagePlus imp2, imp3;
			ImageStack tstack=null, stackNew;
			if (stackSize>1 && doIstack){
				boolean doItAnyway = true;
				if (stackSize>25) {
					YesNoCancelDialog d = new YesNoCancelDialog(IJ.getInstance(),"Auto Local Threshold", "You might run out of memory.\n \nDisplay "+stackSize+" slices?\n \n \'No\' will process without display and\noutput results to the log window.");
					if (!d.yesPressed()){
						//						doIlog=true; //will show in the log window
						doItAnyway=false;
					}
					if (d.cancelPressed())
						return;
				}

				for (int j=1; j<=stackSize; j++){
					imp.setSlice(j);
					ip = imp.getProcessor();
					tstack= new ImageStack(xe,ye);
					for (int k=1; k<ml;k++)
						tstack.addSlice(methods[k], ip.duplicate());
					imp2 = new ImagePlus("Auto Threshold", tstack);
					imp2.updateAndDraw();

					for (int k=1; k<ml;k++){
						imp2.setSlice(k);
						Object[] result = exec(imp2, methods[k], radius, par1, par2, doIwhite );
					}
					//if (doItAnyway){
					CanvasResizer cr= new CanvasResizer();
					stackNew = cr.expandStack(tstack, (xe+2), (ye+18), 1, 1);
					imp3 = new ImagePlus("Auto Threshold", stackNew);
					imp3.updateAndDraw();
					MontageMaker mm= new MontageMaker();
					mm.makeMontage( imp3, 3, 3, 1.0, 1, (ml-1), 1, 0, true); // 3 columns and 3 rows
				}
				imp.setSlice(1);
				//if (doItAnyway)
				IJ.run("Images to Stack", "method=[Copy (center)] title=Montage");
				return;
			}
			else { //single image try all
				tstack= new ImageStack(xe,ye);
				for (int k=1; k<ml;k++)
					tstack.addSlice(methods[k], ip.duplicate());
				imp2 = new ImagePlus("Auto Threshold", tstack);
				imp2.updateAndDraw();

				for (int k=1; k<ml;k++){
					imp2.setSlice(k);
					//IJ.log("analyzing slice with "+methods[k]);
					Object[] result = exec(imp2, methods[k], radius, par1, par2, doIwhite );
				}
				//imp2.setSlice(1);
				CanvasResizer cr= new CanvasResizer();
				stackNew = cr.expandStack(tstack, (xe+2), (ye+18), 1, 1);
				imp3 = new ImagePlus("Auto Threshold", stackNew);
				imp3.updateAndDraw();
				MontageMaker mm= new MontageMaker();
				mm.makeMontage( imp3, 3, 3, 1.0, 1, (ml-1), 1, 0, true);
				return;
			}
		}
		else { // selected a method
			if (stackSize>1 &&  doIstack ) { //whole stack
				//				if (doIstackHistogram) {// one global histogram
				//					Object[] result = exec(imp, myMethod, noWhite, noBlack, doIwhite, doIset, doIlog, doIstackHistogram );
				//				}
				//				else{ // slice by slice
				for (int k=1; k<=stackSize; k++){
					imp.setSlice(k);
					Object[] result = exec(imp, myMethod, radius, par1, par2, doIwhite );
				}
				//				}
				imp.setSlice(1);
			}
			else { //just one slice
				Object[] result = exec(imp, myMethod, radius, par1, par2, doIwhite );
			}
			// 5 - If all went well, show the image:
			// not needed here as the source image is binarised 
		}
	}
	//IJ.showStatus(IJ.d2s((System.currentTimeMillis()-start)/1000.0, 2)+" seconds");


	/** Execute the plugin functionality: duplicate and scale the given image.
	 * @return an Object[] array with the name and the scaled ImagePlus.
	 * Does NOT show the new, image; just returns it. */
	public Object[] exec(ImagePlus imp, String myMethod, int radius,  double par1, double par2, boolean doIwhite ) {

		// 0 - Check validity of parameters
		if (null == imp) return null;
		ImageProcessor ip = imp.getProcessor();
		int xe = ip.getWidth();
		int ye = ip.getHeight();

		//int [] data = (ip.getHistogram());

		IJ.showStatus("Thresholding...");
		long startTime = System.currentTimeMillis();
		//1 Do it
		if (imp.getStackSize()==1){
			ip.snapshot();
			Undo.setup(Undo.FILTER, imp);
		}
		// Apply the selected algorithm
		if(myMethod.equals("Bernsen")){
			Bernsen(imp,  radius, par1, par2, doIwhite);
		}
		else if(myMethod.equals("BernsenOtsu")){
			BernsenOtsu(imp, radius, par1, par2, doIwhite);
		}
		else if(myMethod.equals("std_dev")){
			std(imp, radius, par1, par2, doIwhite);
		}
		else if(myMethod.equals("stdwithBernsen")){
			stdwithBernsen(imp, radius, par1, par2, doIwhite);
		}
		else if(myMethod.equals("stdwithBernsenKim")){
			stdwithBernsenKim(imp, radius, par1, par2, doIwhite);
		}
		//IJ.showProgress((double)(255-i)/255);
		imp.updateAndDraw();
		imp.getProcessor().setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE);
		// 2 - Return the threshold and the image
		IJ.showStatus("\nDone " + (System.currentTimeMillis() - startTime) / 1000.0);
		return new Object[] {imp};
	}


	private void std(ImagePlus imp, int radius, double par1, double par2, boolean doIwhite) {
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



	void stdwithBernsen(ImagePlus imp, int radius, double par1, double par2, boolean doIwhite) {
		int w=imp.getWidth();
		int h=imp.getHeight();
		int position;
		int radiusx2=radius * 2;
		float object;
		float backg;

		if (doIwhite){
			object =  (float) 0xff;
			backg =   (float) 0;
		}
		else {
			object =  (float) 0;
			backg =  (float) 0xff;
		}
		
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
				float local_contrast=localmax-localmin;
				float midgray=(localmax+localmin)/2;
				float k=(float)5;
				
				//standard deviation will become local threshold
				if ( localstd < k )
					copy[position] = ( midgray >= 128.0 ) ? object :  backg;  //Low contrast region
				else
					copy[position] = (pixels[position] > midgray ) ? object : backg;
			}
		}
		for (position=0; position<w*h; position++) pixels[position]=copy[position];
	}


	void BernsenOtsu(ImagePlus imp, int radius, double par1, double par2, boolean doIwhite) {
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

		if (doIwhite){
			object =  (byte) 0xff;
			backg =   (byte) 0;
		}
		else {
			object =  (byte) 0;
			backg =  (byte) 0xff;
		}

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

				if (local_contrast < 15)
					pixelsOut[position] = ( midgray >= 128 ) ? object :  backg;
				else	
					pixelsOut[position] = ((int) (pixels[position]&0xff)>kStar) ? object : backg;
			}
		}
		for (position=0; position<w*h; position++) pixels[position]=pixelsOut[position]; //update with thresholded pixels

	}

	void stdwithBernsenKim(ImagePlus imp, int radius, double par1, double par2, boolean doIwhite) {
		int w=imp.getWidth();
		int h=imp.getHeight();
		int radiusx2=radius * 2;

		
		ImagePlus curr=imp;
		ImageConverter convrt=new ImageConverter(curr);
		convrt.convertToGray32();

		ImageProcessor ip=curr.getProcessor();
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
				if (localstd < 5) {
					numlowstds++;
					sumlowstds=sumlowstds+localstd;
				}
				
			}
		}
		double meanlowstds=sumlowstds/(double)numlowstds;
		double k=18;
		double local_threshold=k*meanlowstds;
		Bernsen32(curr, radius, local_threshold, par2, doIwhite);
	}
	
	void Bernsen32(ImagePlus imp, int radius,  double par1, double par2, boolean doIwhite ) {
		int w=imp.getWidth();
		int h=imp.getHeight();
		int position;
		int radiusx2=radius * 2;
		float object;
		float backg;
		int contrast_threshold=15;
		
		if (par1!=0) {
			IJ.log("Bernsen: changed contrast_threshold from :"+ contrast_threshold + "  to:" + par1);
			contrast_threshold= (int) par1;
		}
		
		if (doIwhite){
			object =  (float) 0xff;
			backg =   (float) 0;
		}
		else {
			object =  (float) 0;
			backg =  (float) 0xff;
		}
		
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


				float localmax=(float) ip.getStats().max;
				float localmin=(float) ip.getStats().min;
				float local_contrast=localmax-localmin;
				float midgray=(localmax+localmin)/2;
				
				if ( local_contrast < contrast_threshold )
					copy[position] = ( midgray > 128 ) ? object :  backg;  //Low contrast region
				else
					copy[position] = (pixels[position] >= midgray ) ? object : backg;
			}
		}
		for (position=0; position<w*h; position++) pixels[position]=copy[position];
	}
	
	void Bernsen(ImagePlus imp, int radius,  double par1, double par2, boolean doIwhite ) {
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
		int contrast_threshold=15;
		int local_contrast;
		int mid_gray;
		byte object;
		byte backg;
		int temp;

		if (par1!=0) {
			IJ.log("Bernsen: changed contrast_threshold from :"+ contrast_threshold + "  to:" + par1);
			contrast_threshold= (int) par1;
		}

		if (doIwhite){
			object =  (byte) 0xff;
			backg =   (byte) 0;
		}
		else {
			object =  (byte) 0;
			backg =  (byte) 0xff;
		}

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
				pixels[i] = ( mid_gray >= 128 ) ? object :  backg;  //Low contrast region
			else
				pixels[i] = (temp >= mid_gray ) ? object : backg;
		}    
		//imp.updateAndDraw();
		return;
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
