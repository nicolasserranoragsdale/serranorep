import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.StackWindow;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.util.ArrayUtil;

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
	private ImageStack stack;
	private int c;
	private int bitDepth;
 
	public int setup(String arg, ImagePlus imp) {
		imRef = imp;
		stack=imp.getImageStack();
		
		bitDepth=imp.getBitDepth();
		
		if (arg.equals("about")) {
			showAbout();
			return DONE;
		}

		getParams();

		return DOES_ALL;
	}

	private void getParams() {

		// set defaults
		if (bitDepth == 16) {
			baseThreshold = 40000;
			k = 500;
		}
		else {
			baseThreshold = 128;
			k = 5;
		}
		
		int diameter = 10;
		c=18;

		GenericDialog gd = new GenericDialog(
				"3D threshold");

		gd.addNumericField("Base threshold", baseThreshold, 0);
		gd.addNumericField("Mask diameter (pixels)", diameter, 0);
		gd.addNumericField("low std_dev threshold", k, 1);
		gd.addNumericField("multiplication factor for std_dev contrast threshold", c, 0);

		String [] methods={"std_dev_plot","BernsenKim3Dellipsoid","BernsenKim3Dsphere"};
		gd.addChoice("Method", methods, methods[1]);

		gd.showDialog();

		if (gd.wasCanceled()) {
			if (imRef != null)
				imRef.unlock();
			noGo = true;
		}

		baseThreshold = (int) gd.getNextNumber();
		radius = ((int) gd.getNextNumber()) / 2;
		k = gd.getNextNumber() ;
		c=(int)gd.getNextNumber() ;
		myMethod=gd.getNextChoice ();
	}

	public void run(ImageProcessor ip) {

		if (noGo)
			return;
		if(myMethod.equals("std_dev_plot")){
			std(ip);
		}
		else if(myMethod.equals("BernsenKim3Dsphere")){
			BernsenKimsphere(ip);
		}
		else if(myMethod.equals("BernsenKim3Dellipsoid")){
			BernsenKimellipsoid(ip);
		}
	}




	private void std(ImageProcessor ip) {
		width = ip.getWidth();
		height = ip.getHeight();
		depth = imRef.getStackSize();
		

		int radx=radius;
		int rady=radius;
		int radz=radius;
		int zmin=0;
		int zmax=depth;
		float local_std_dev=0;
		
		int[] ker = createKernelEllipsoid(radx, rady, radz);
        int nb = 0;
        for (int i=0; i<ker.length; i++)
            nb += ker[i];
        
        
        // create variable to store the new, binary image
        float[][] imageCopy = new float[depth][width * height];
        
        for (int z=zmin; z<zmax; z++) {
            if (zmin==0) IJ.showProgress(z+1, zmax);
            for (int y=0; y<height; y++) {
                for (int x=0; x<width; x++) {
                    ArrayUtil tab = getNeighborhood(ker, nb, x, y, z, radx, rady, radz);
                    local_std_dev=(float) Math.sqrt(tab.getVariance());
                    
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
	
	private void BernsenKimsphere(ImageProcessor ip) {
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
		
	}
	
	
	
	
	private void BernsenKimellipsoid(ImageProcessor ip) {
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

		int radx=radius;
		int rady=radius;
		int radz=radius;
		int zmin=0;
		int zmax=depth;
		
		int[][] local_contrast=new int[depth][width * height];
		float[][] mid_gray=new float[depth][width*height];
		
		int[] ker = createKernelEllipsoid(radx, rady, radz);
        int nb = 0;
        for (int i=0; i<ker.length; i++)
            nb += ker[i];
        
        for (int z=zmin; z<zmax; z++) {
            if (zmin==0) IJ.showProgress(z+1, zmax);
            for (int y=0; y<height; y++) {
                for (int x=0; x<width; x++) {
                    ArrayUtil tab = getNeighborhood(ker, nb, x, y, z, radx, rady, radz);
                    
                    int localMax=(int) tab.getMaximum();
                    int localMin=(int) tab.getMinimum();
                    local_contrast[z ][x + y * width]=localMax-localMin;
					mid_gray[z][x + y * width]=(float)(localMax+localMin)/(float)2.0;
					local_std_dev=Math.sqrt(tab.getVariance());
					if (local_std_dev < k) {
						sumlowstds+=local_std_dev;
						countlowstds++;
					}
                    
                }
            }
        }
        meanlowstds=sumlowstds/(double)countlowstds;
		
	
		double contrastThreshold=c*meanlowstds;
		//double contrastThreshold=15;
		IJ.log("contrast threshold="+ contrastThreshold);
		
		for (int z=zmin; z<zmax; z++) {
            IJ.showProgress(z+1, zmax);
            for (int y=0; y<height; y++) {
                for (int x=0; x<width; x++) {
                	
                	
                	int value = (int)stack.getVoxel(x, y, z);
                	
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

		ImagePlus newImage = new ImagePlus("3DThresholdEllipse", newStack);
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
	
	
	
	/**
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

	void showAbout() {
		IJ.showMessage("About Adaptive 3D Threshold..",
				"This plugin thresholds a stack according to the threshold "
						+ "T = (1-w)*base - w*avg(radius^3 neighbour-base)");
	}

}
