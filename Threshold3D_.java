import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ProgressBar;
import ij.gui.StackWindow;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class Threshold3D_ implements PlugInFilter {
	ImagePlus imRef;
	private boolean noGo;
	private int baseThreshold;
	private int radius;
	private double k;
	private int width;
	private int height;
	private int depth;

	private String myMethod;

	public int setup(String arg, ImagePlus imp) {
		imRef = imp;

		if (arg.equals("about")) {
			showAbout();
			return DONE;
		}

		getParams();

		return DOES_8G;
	}

	private void getParams() {

		// set defaults
		baseThreshold = 128;
		k = 5;
		int diameter = 3;

		GenericDialog gd = new GenericDialog(
				"3D threshold");

		gd.addNumericField("Base threshold", baseThreshold, 0);
		gd.addNumericField("Mask diameter (pixels)", diameter, 0);
		gd.addNumericField("std_dev threshold", k, 1);

		String [] methods={"std_dev_plot","BernsenKim3D"};
		gd.addChoice("Method", methods, methods[0]);

		gd.showDialog();

		if (gd.wasCanceled()) {
			if (imRef != null)
				imRef.unlock();
			noGo = true;
		}

		baseThreshold = (int) gd.getNextNumber();
		radius = ((int) gd.getNextNumber()) / 2;
		k = gd.getNextNumber() ;
		myMethod=gd.getNextChoice ();
	}

	public void run(ImageProcessor ip) {

		if (noGo)
			return;
		if(myMethod.equals("std_dev_plot")){
			std(ip);
		}
		else if(myMethod.equals("BernsenKim3D")){
			BernsenKim(ip);
		}
	}



	private void std(ImageProcessor ip) {

		width = ip.getWidth();
		height = ip.getHeight();
		depth = imRef.getStackSize();

		// create variable to store the new, binary image
		float[][] imageCopy = new float[depth][width * height];

		// do the thresholding
		int localValue;
		long sum;
		long count;
		double localAverage;

		float local_std_dev=0;
		for (int z = 0; z < depth; z++) {
			IJ.showProgress(z , depth);
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {

					sum = 0;
					count = 0;
					localAverage = 0;
					for (int i = -radius; i < radius; i++) {
						for (int j = -radius; j < radius; j++) {
							for (int k1 = -radius; k1 < radius; k1++) {
								localValue = safeGet(z + k1 + 1, x + i, y + j);

								if (localValue >= 0) { // if not outside mask
									sum += localValue;
									count++;
								}
							}
						}
					}


					localAverage = sum / count;

					//traverse again to calculate std_dev using local average
					double sumOfSquares=0;
					local_std_dev=0;
					for (int i = -radius; i < radius; i++) {
						for (int j = -radius; j < radius; j++) {
							for (int k1 = -radius; k1 < radius; k1++) {
								localValue = safeGet(z + k1 + 1, x + i, y + j);

								if (localValue >= 0)  // if not outside mask
									sumOfSquares += Math.pow((localValue-localAverage),2);

							}
						}
					}



					local_std_dev=(float)Math.sqrt(sumOfSquares/(count-1));

					imageCopy[z][x + y * width] = local_std_dev;
				}
			}

		}
		ImageStack newStack = new ImageStack(width, height);
		
		
		for (int i = 0; i < depth; i++) {
			float[] newPixels = imageCopy[i];
			newStack.addSlice("Slice " + i, newPixels);
		}


		ImagePlus newImage = new ImagePlus("3DThreshold", newStack);
		new StackWindow(newImage);

	}


	private void BernsenKim(ImageProcessor ip) {
		
		width = ip.getWidth();
		height = ip.getHeight();
		depth = imRef.getStackSize();

		// create variable to store the new, binary image
		byte[][] imageCopy = new byte[depth][width * height];

		// do the thresholding
		int value;
		int localValue;
		long sum;
		long count;
		double contrastThreshold;
		double localAverage;

		double local_std_dev=0;
		double meanlowstds=0;
		double sumlowstds=0;
		int countlowstds=0;
		
		int[][] local_contrast=new int[depth][width * height];
		float[][] mid_gray=new float[depth][width*height];
		for (int z = 0; z < depth; z++) {
			IJ.showProgress(z + 1, depth/2);
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {

					sum = 0;
					count = 0;
					localAverage = 0;
					for (int i = -radius; i < radius; i++) {
						for (int j = -radius; j < radius; j++) {
							for (int k1 = -radius; k1 < radius; k1++) {
								localValue = safeGet(z + k1 + 1, x + i, y + j);

								if (localValue >= 0) { // if not outside image stack
									sum += localValue;
									count++;
								}
							}
						}
					}


					localAverage = sum / count;

					//traverse again to calculate std_dev using local average
					double sumOfSquares=0;
					local_std_dev=0;
					int localMax=0;
					int localMin=255;
					for (int i = -radius; i < radius; i++) {
						for (int j = -radius; j < radius; j++) {
							for (int k1 = -radius; k1 < radius; k1++) {
								localValue = safeGet(z + k1 + 1, x + i, y + j);
								
								if (localValue>localMax)
									localMax=localValue;
								if (localValue<localMin)
									localMin=localValue;
								
								if (localValue >= 0)  // if not outside mask
									sumOfSquares += Math.pow((localValue-localAverage),2);

							}
						}
					}
					
					local_contrast[z ][x + y * width]=localMax-localMin;
					mid_gray[z][x + y * width]=(float)(localMax+localMin)/2;
					
					local_std_dev=Math.sqrt(sumOfSquares/(count-1));
					if (local_std_dev < k) {
						sumlowstds+=local_std_dev;
						countlowstds++;
					}
					
					/*localThreshold = (int) (baseThreshold - k
							* localAverage + 0.5);

					if (value >= localThreshold) {
						imageCopy[z][x + y * width] = (byte) 255;
					} else {
						imageCopy[z][x + y * width] = (byte) 0;
					}*/
				}
			}

		}
		
		
		
		meanlowstds=sumlowstds/countlowstds;
		
		double c=18;
		contrastThreshold=c*meanlowstds; 
		//use low standard deviations to determine contrast threshold in bernsten
		IJ.showMessage("madeit");
		
		for (int z = 0; z < depth; z++) {
			IJ.showProgress(depth/2+z, depth - 1);
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					
					value = 0xff & ((byte[]) imRef.getStack().getPixels(z + 1))[x+ y * width];
					
					if ( local_contrast[z][x + y * width] < contrastThreshold )
						imageCopy[z][x + y * width] = 
							(mid_gray[z][x + y * width] >= baseThreshold ) ? (byte) 255 :  (byte) 0;  //Low contrast region
					else
						imageCopy[z][x + y * width] = 
							(value >= mid_gray[z][x + y * width] ) ? (byte) 255 :  (byte) 0;
	
				}
			}
		}
		
		
		
		

		ImageStack newStack = new ImageStack(width, height);

		for (int i = 0; i < depth; i++) {
			byte[] newPixels = imageCopy[i];
			newStack.addSlice("Slice " + i, newPixels);
		}

		IJ.showProgress(1, 1); // set to finished

		ImagePlus newImage = new ImagePlus("3DThreshold", newStack);
		new StackWindow(newImage);



	}



	private int safeGet(int z, int x, int y) {

		// Gets the value from the image, or if outside image return -1

		int retval;

		try {
			retval = 0xff & ((byte[]) imRef.getStack().getPixels(z + 1))[x + y
			                                                             * width];
		} catch (Exception e) {
			retval = -1;

		}

		return retval;

	}

	void showAbout() {
		IJ.showMessage("About Adaptive 3D Threshold..",
				"This plugin thresholds a stack according to the threshold "
						+ "T = (1-w)*base - w*avg(radius^3 neighbour-base)");
	}

}
