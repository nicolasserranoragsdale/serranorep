import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.StackWindow;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.util.ArrayUtil;
//Class contains a few 3D adaptive thresholding methods for 8 and 16 bit greyscale images, based on
//Bernsen.
public class Threshold3D_ implements PlugInFilter {
	ImagePlus imRef;
	private boolean noGo;
	private int baseThreshold;
	private int radius;
	private int contrastThreshold;
	private int width;
	private int height;
	private int depth;

	private String myMethod;
	private ImageStack stack;
	private int c;
	private double c2;
	private int bitDepth;

	public int setup(String arg, ImagePlus imp) {
		imRef = imp;
		stack=imp.getImageStack();

		bitDepth=imp.getBitDepth();

		getParams();

		return DOES_ALL;
	}

	private void getParams() {

		// set defaults depending on bit depth
		if (bitDepth == 16) {
			baseThreshold = 40000;
			contrastThreshold=2500;
		}
		else if (bitDepth == 8) {
			baseThreshold = 128;
			contrastThreshold=15;
		}

		radius = 5;
		//multiplication factors
		c=18;
		c2=6;

		
		//Create generic dialog to get user input values
		GenericDialog gd = new GenericDialog(
				"3D threshold");

		gd.addNumericField("Base threshold", baseThreshold, 0);
		gd.addNumericField("Mask radius (pixels)", radius, 0);
		gd.addNumericField("Contrast Threshold", contrastThreshold, 0);
		gd.addNumericField("multiplication factor (BernsenKim)", c, 0);
		gd.addNumericField("multiplication factor (BersenContrast)", c2, 0);

		String [] methods={"BernsenKim3Dellipsoid",
				"BernsenContrast3D"};
		
		gd.addChoice("Method", methods, methods[0]);

		gd.showDialog();
		
		if (gd.wasCanceled()) {
			if (imRef != null)
				imRef.unlock();
			noGo = true;
		}

		
		//set values based on user info
		baseThreshold = (int) gd.getNextNumber();
		radius = ((int) gd.getNextNumber()) ;
		contrastThreshold=(int)gd.getNextNumber();
		c=(int)gd.getNextNumber() ;
		c2=gd.getNextNumber();
		myMethod=gd.getNextChoice ();
	}

	public void run(ImageProcessor ip) {
		if (noGo)
			return;
		
		long start = System.currentTimeMillis();
		
		//choose desired method 
		if(myMethod.equals("BernsenKim3Dellipsoid")){
			BernsenKimellipsoid(ip);
		}
		else if(myMethod.equals("BernsenContrast3D")){
			BernsenContrast3D(ip);
		}
		
		long end = System.currentTimeMillis();
		//will keep track of how long method takes
		IJ.log("time(seconds): "+(double)((end-start)/(double)1000));
	}


	//Method performs a Bernsen-style thresholding algorithm on image by looking at spherical region 
	//around each pixel. First travels through stack to store local contrast and midgray values in a 
	//large array and to calculate a fair local threshold value based on std_dev of low contrast regions, 
	//then sets each pixel to white or black accordingly.
	//Note: Kernel not actually ellipsoid. It is a sphere, but can easily be modified to be ellipsoid.
	private void BernsenKimellipsoid(ImageProcessor ip) {
		//get dimensions
		width = ip.getWidth();
		height = ip.getHeight();
		depth = imRef.getStackSize();
		stack=imRef.getImageStack();

		double local_std_dev=0;
		double meanlowstds=0;
		double sumlowstds=0;
		int countlowstds=0;
		
		// create variable to store the new, binary image
		byte[][] imageCopy = new byte[depth][width * height];
		
		//made into sphere for now
		int radx=radius;
		int rady=radius;
		int radz=radius;
		int zmin=0;
		int zmax=depth;

		//intializing large arrays to characterize spherical regions around each pixel
		int[][] local_contrast=new int[depth][width * height];
		float[][] mid_gray=new float[depth][width*height];

		//calls method to store kernel as a boolean array
		int[] ker = createKernelEllipsoid(radx, rady, radz);
		int nb = 0;
		for (int i=0; i<ker.length; i++)
			nb += ker[i];
		
		//travel through all pixels in stack
		for (int z=zmin; z<zmax; z++) {
			if (zmin==0) IJ.showProgress(z+1, zmax);
			for (int y=0; y<height; y++) {
				for (int x=0; x<width; x++) {
					//gets surrounding pixels
					ArrayUtil tab = getNeighborhood(ker, nb, x, y, z, radx, rady, radz);

					int localMax=(int) tab.getMaximum();
					int localMin=(int) tab.getMinimum();
					//store local contrast and midgray values in array using local max and min.
					local_contrast[z ][x + y * width]=localMax-localMin;
					mid_gray[z][x + y * width]=(float)(localMax+localMin)/(float)2.0;
					
					//if local contrast is less than the user-defined threshold and region is within
					//bounds, consider region's std_dev
					if (((localMax-localMin) < contrastThreshold) && localMax!=0 && localMin!=0) {
						local_std_dev=Math.sqrt(tab.getVariance());
						sumlowstds+=local_std_dev;
						countlowstds++;
					}

				}
			}
		}
		//calculate the mean standard deviation of low contrast regions
		meanlowstds=sumlowstds/(double)countlowstds;

		//mean std_dev times some constant defined by user will become new contrast threshold
		double contrastThreshold=c*meanlowstds;
		IJ.log("contrast threshold="+ contrastThreshold);

		//travel through stack again to set pixels to white or black
		for (int z=zmin; z<zmax; z++) {
			//variables to keep track of porosity
			int porecount=0;
			int whitecount=0;
			IJ.showProgress(z+1, zmax);
			for (int y=0; y<height; y++) {
				for (int x=0; x<width; x++) {

					//value of current pixel
					int value = (int)stack.getVoxel(x, y, z);
					
					//BERNSEN
					//First evaluates if region is high or low contrast, then changes pixel 
					if ( local_contrast[z][x + y * width] < contrastThreshold ) { //Low contrast region
						//make sure pixel is within valid bounds
						if (mid_gray[z][x + y * width] >= baseThreshold && value!=0) {
							whitecount++;
							imageCopy[z][x + y * width] =(byte) 255; //set to white
						}
						if (mid_gray[z][x + y * width] < baseThreshold && value!=0) {
							porecount++;
							imageCopy[z][x + y * width] =(byte) 0; //set to black
						}
						//imageCopy[z][x + y * width] = 
						//	(mid_gray[z][x + y * width] >= baseThreshold ) ? (byte) 255 :  (byte) 0;  //Low contrast region
					}
					else { //High contrast region
						if (value >= mid_gray[z][x + y * width] && value!=0) {
							whitecount++;
							imageCopy[z][x + y * width] =(byte) 255;
						}
						if (value < mid_gray[z][x + y * width] && value!=0) {
							porecount++;
							imageCopy[z][x + y * width] =(byte) 0;
						}
						//imageCopy[z][x + y * width] = 
						//	(value >= mid_gray[z][x + y * width] ) ? (byte) 255 :  (byte) 0;
					}
				}
			}
			//porosity level for slice
			IJ.log(""+(float)(porecount)*100/(float)(whitecount+porecount));
		}
		//displays results in new stack
		ImageStack newStack = new ImageStack(width, height);

		for (int i = 0; i < depth; i++) {
			byte[] newPixels = imageCopy[i];
			newStack.addSlice("Slice " + i, newPixels);
		}

		IJ.showProgress(1, 1); // set to finished

		ImagePlus newImage = new ImagePlus("3DBernsenKim", newStack);
		new StackWindow(newImage);

	}


	//Very similar to BernsenKim. It performs a Bernsen-style thresholding algorithm on image by looking
	//at spherical region around each pixel. 
	//However, when traveling through stack to store values, it calculates the mean and std_dev
	//of the localcontrast value, instead of the pixels within the region themselves.
	private void BernsenContrast3D(ImageProcessor ip) {
		//get dimensions
		width = ip.getWidth();
		height = ip.getHeight();
		depth = imRef.getStackSize();
		stack=imRef.getImageStack();

		
		// create variable to store the new, binary image
		byte[][] imageCopy = new byte[depth][width * height];

		//made into sphere for now
		int radx=radius;
		int rady=radius;
		int radz=radius;
		int zmin=0;
		int zmax=depth;

		//intializing large arrays to characterize spherical regions around each pixel
		int[][] local_contrast=new int[depth][width * height];
		float[][] mid_gray=new float[depth][width*height];

		//calls method to store kernel as a boolean array
		int[] ker = createKernelEllipsoid(radx, rady, radz);
		int nb = 0;
		for (int i=0; i<ker.length; i++)
			nb += ker[i];
		
		
		
		int sum_lowcon=0;
		int num_lowcon=0;
		
		//travel through all pixels in stack
		for (int z=zmin; z<zmax; z++) {
			if (zmin==0) IJ.showProgress(z+1, zmax);
			for (int y=0; y<height; y++) {
				for (int x=0; x<width; x++) {
					//gets surrounding pixels
					ArrayUtil tab = getNeighborhood(ker, nb, x, y, z, radx, rady, radz);

					int localMax=(int) tab.getMaximum();
					int localMin=(int) tab.getMinimum();
					//store local contrast and midgray values in array using local max and min.
					local_contrast[z ][x + y * width]=localMax-localMin;
					mid_gray[z][x + y * width]=(float)(localMax+localMin)/(float)2.0;
			
					//if local contrast is less than the user-defined threshold and region is within
					//bounds, consider region's local contrast
					if (((localMax-localMin) < contrastThreshold) && (localMax-localMin!=0)) {
						sum_lowcon+=(localMax-localMin);
						num_lowcon++;
					}

				}
			}
		}
		//calculate average local contrast for the low contrast regions
		double mean_lowcon=(double)sum_lowcon/(double)num_lowcon;
		
		//travel through array to calculate std_dev of the localcontrasts in low contrast regions
		double sumsquares_lowcon=0;
		for (int z=zmin; z<zmax; z++) {
			if (zmin==0) IJ.showProgress(z+1, zmax);
			for (int y=0; y<height; y++) {
				for (int x=0; x<width; x++) {
					 
					//check to see if low contrast and then if it's a valid region
					if (local_contrast[z ][x + y * width] < contrastThreshold 
							&& local_contrast[z ][x + y * width] !=0) {
						sumsquares_lowcon+=(local_contrast[z ][x + y * width]-mean_lowcon)
								*(local_contrast[z ][x + y * width]-mean_lowcon);
					}

				}
			}
		}
		double std_lowcon=Math.sqrt(sumsquares_lowcon/((double)num_lowcon-1));
		
		//Mean localcontrast plus the std_dev of localcontrasts times a constant becomes new contrast
		//threshold.
		double contrastThreshold=mean_lowcon+c2*std_lowcon;
		//double contrastThreshold=15;
		//IJ.log("contrast threshold="+ contrastThreshold);


		//travel through stack again to set pixels to white or black
		for (int z=zmin; z<zmax; z++) {
			//variables to keep track of porosity
			int porecount=0;
			int whitecount=0;
			IJ.showProgress(z+1, zmax);
			for (int y=0; y<height; y++) {
				for (int x=0; x<width; x++) {


					int value = (int)stack.getVoxel(x, y, z);

					//BERNSEN
					//First evaluates if region is high or low contrast, then changes pixel 
					if ( local_contrast[z][x + y * width] < contrastThreshold ) { //Low contrast region
						if (mid_gray[z][x + y * width] >= baseThreshold && value!=0) {
							whitecount++;
							imageCopy[z][x + y * width] =(byte) 255; //set to white
						}
						if (mid_gray[z][x + y * width] < baseThreshold && value!=0) {
							porecount++;
							imageCopy[z][x + y * width] =(byte) 0;  //set to black
						}
						//imageCopy[z][x + y * width] = 
						//	(mid_gray[z][x + y * width] >= baseThreshold ) ? (byte) 255 :  (byte) 0;  //Low contrast region
					}
					else {	      //High contrast region
						if (value >= mid_gray[z][x + y * width] && value!=0) { 
							whitecount++;
							imageCopy[z][x + y * width] =(byte) 255;
						}
						if (value < mid_gray[z][x + y * width] && value!=0) {
							porecount++;
							imageCopy[z][x + y * width] =(byte) 0;
						}
						//imageCopy[z][x + y * width] = 
						//	(value >= mid_gray[z][x + y * width] ) ? (byte) 255 :  (byte) 0;
					}
				}
			}
			//porosity level for slice
			IJ.log(""+(float)(porecount)*100/(float)(whitecount+porecount));
		}
		ImageStack newStack = new ImageStack(width, height);

		for (int i = 0; i < depth; i++) {
			byte[] newPixels = imageCopy[i];
			newStack.addSlice("Slice " + i, newPixels);
		}

		IJ.showProgress(1, 1); // set to finished

		ImagePlus newImage = new ImagePlus("3DBernsenContrast", newStack);
		new StackWindow(newImage);

	}
	
	
	
	
	
	




	/**
	 * Method increases efficiency by creating an array of 1 and 0's that will show us which pixels are
	 * within the ellipsoid. Created by Thomas Boudier for adaptive filtering.
	 * 
	 * Thomas Boudier Create a kernel neighorhood as an ellipsoid
	 *
	 * @param radx Radius x of the ellipsoid
	 * @param rady Radius x of the ellipsoid
	 * @param radz Radius x of the ellipsoid
	 * @return The kernel as an array
	 */
	private int[] createKernelEllipsoid(float radx, float rady, float radz) {
		int vx = (int) Math.ceil(radx);
		int vy = (int) Math.ceil(rady);
		int vz = (int) Math.ceil(radz);
		int[] ker = new int[(2 * vx + 1) * (2 * vy + 1) * (2 * vz + 1)];
		double dist;

		double rx2 = radx * radx;
		double ry2 = rady * rady;
		double rz2 = radz * radz;

		if (rx2 != 0) {
			rx2 = 1.0 / rx2;
		} else {
			rx2 = 0;
		}
		if (ry2 != 0) {
			ry2 = 1.0 / ry2;
		} else {
			ry2 = 0;
		}
		if (rz2 != 0) {
			rz2 = 1.0 / rz2;
		} else {
			rz2 = 0;
		}

		int idx = 0;
		for (int k = -vz; k <= vz; k++) {
			for (int j = -vy; j <= vy; j++) {
				for (int i = -vx; i <= vx; i++) {
					dist = ((double) (i * i)) * rx2 + ((double) (j * j)) * ry2 + ((double) (k * k)) * rz2;
					if (dist <= 1.0) {
						ker[idx] = 1;
					} else {
						ker[idx] = 0;
					}
					idx++;
				}
			}
		}

		return ker;
	}


	//Also created by Thomas Boudier for adaptive filtering. Returns voxels within ellipsoid as an 
	//ArrayUtil by checking kernel ellipsoid array.
	private ArrayUtil getNeighborhood(int[] ker, int nbval, int x, int y, int z, float radx, float rady, float radz) {
		ArrayUtil pix = new ArrayUtil(nbval);
		int vx = (int) Math.ceil(radx);
		int vy = (int) Math.ceil(rady);
		int vz = (int) Math.ceil(radz);
		int index = 0;
		int c = 0;
		int sizex = stack.getWidth();
		int sizey = stack.getHeight();
		int sizez = stack.getSize();
		for (int k = z - vz; k <= z + vz; k++) {
			for (int j = y - vy; j <= y + vy; j++) {
				for (int i = x - vx; i <= x + vx; i++) {
					if (ker[c]>0 && i>=0 && j>=0 && k>=0 && i<sizex && j<sizey && k<sizez) {
						pix.putValue(index, (float)stack.getVoxel(i, j, k));
						index++;
					}
					c++;
				}
			}
		}
		pix.setSize(index);
		return pix;
	}

	
	
	
	//previous method inefficient. No longer in use.
	
	
	
/*	 * $Id: Adapative3DThreshold_.java,v 1.22 2005/07/05 12:40:32 perchrh Exp $
	 * 
	 * Created by Per Christian Henden. Free software in the public domain.*/
	
		/*private void BernsenKimsphere(ImageProcessor ip) {
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
				IJ.showProgress(z + 1, depth*2);
				for (int y = 0; y < height; y++) {
					for (int x = 0; x < width; x++) {

						sum = 0;
						count = 0;
						localAverage = 0;
						for (int i = -radius; i < radius; i++) {
							for (int j = -radius; j < radius; j++) {
								for (int k1 = -radius; k1 < radius; k1++) {
									if(radius*radius>=i*i+j*j+k1*k1) { //check if in sphere
										localValue = safeGet(z + k1 + 1, x + i, y + j);

										if (localValue >= 0) { // if not outside image stack
											sum += localValue;
											count++;
										}
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
									if(radius*radius>=i*i+j*j+k1*k1) { //check if in sphere
										localValue = safeGet(z + k1 + 1, x + i, y + j);

										if (localValue >= 0) { // if not outside mask

											if (localValue>localMax)
												localMax=localValue;
											if (localValue<localMin)
												localMin=localValue;

											sumOfSquares += Math.pow((localValue-localAverage),2);
										}
									}

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

						localThreshold = (int) (baseThreshold - k
						 * localAverage + 0.5);

						if (value >= localThreshold) {
							imageCopy[z][x + y * width] = (byte) 255;
						} else {
							imageCopy[z][x + y * width] = (byte) 0;
						}
					}
				}

			}



			meanlowstds=sumlowstds/countlowstds;


			contrastThreshold=c*meanlowstds; 
			//use low standard deviations to determine contrast threshold in bernsten

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

			IJ.showProgress(1, 1); 

			ImagePlus newImage = new ImagePlus("3DThreshold", newStack);
			new StackWindow(newImage);

		}*/
	
	
	
	/*private int safeGet(int z, int x, int y) {

	// Gets the value from the image, or if outside image return -1

	int retval;

	try {
		retval = 0xff & ((byte[]) imRef.getStack().getPixels(z + 1))[x + y
		                                                             * width];
	} catch (Exception e) {
		retval = -1;

	}

	return retval;

}*/

}
