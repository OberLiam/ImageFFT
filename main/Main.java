package main;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Random;

import javax.imageio.ImageIO;


/**
 * @author Liam Axon
 * Overall notes:
 *     must set w, h, inFile, outFile to use
 *     if using FFT, must have w, h powers of 2
 *     handles RGB values at face value (i.e. no conversion to brightness levels)
 *     will ignore all alpha values
 */
public class Main {
	
	
	//NOTE: w, h used throughout. Not to be changed
	public static int W = 1024; //width
	public static int H = 512; //height
	public static int M = 2048; //common multiple of H and W used for FFT
	public static String INFILE = "imgs/chicago.png"; // the file to read from (include ".png")
	public static String OUTFILE = "imgs/test/chicago"; // the file to write to (no ".png")
	
	public static Random random = new Random(); //utility random instance used everywhere
	
	public static void main(String[] args) {
		System.out.println("Starting.");
		
		// Setting up expsM as the roots of unity
		expsM = new Complex[M];
		for(int i=0;i<M;i++) {
			expsM[i] = Complex.makeRT(1, (2*Math.PI)*i/(double)(M) );
		}
		
		// Setting up the random variable generator
		random = new Random();
		
		// Now the main part of the program:
		try{
			System.out.println("Reading in picture.");
			// Reading in a picture, the converting it to wave-form:
			BufferedImage bimg = ImageIO.read(new File(INFILE));
			Complex[][][] pxArr = getComplexPxlArr(bimg);
			double targetBrightness = getBrightness(pxArr);
			
			System.out.println("Performing Fourier Transform...");
			long temp1 = System.currentTimeMillis();
			Complex[][][] wvArr = FTAll(pxArr, true);
			System.out.println("Finished Fourier Transform.");
			System.out.println("(took "+(System.currentTimeMillis()-temp1)/1000.0+" seconds)");

//			//Writing the picture (a pure copy):
//			BufferedImage outBImg = makeBImg(pxArr);
//			ImageIO.write(outBImg, "png", new File(OUTFILE+"_copy.png"));
//			System.out.println("Finished printing 'copy'");
			
//			//Testing that the Fourier transform does invert itself:
//			Complex[][][] pxArrTest = FFTAll(wvArr, false);
//			//Writing the picture (test).
//			BufferedImage testBImg = makeBImg(pxArrTest);
//			ImageIO.write(testBImg, "png", new File(OUTFILE+"_test.png"));
//			System.out.println("Finished printing 'test'");
			
			
//			//Here's the low part:
//			Filterer filterLow = (W, H, i, j, clr) -> i+j < (W+H)/2;
//			Complex[][][] pxArrLow = FFTAll(applyFilter(wvArr, filterLow), false);
//			ImageIO.write( makeBImg(pxArrLow), "png", new File(OUTFILE+"_low.png"));
			
//			//Here's the high part
//			Filterer filterHigh = (W, H, i, j, clr) -> i+j > 0;
//			Complex[][][] pxArrHigh = FFTAll(applyFilter(wvArr, filterHigh), false);
//			makeBright(pxArrHigh, targetBrightness);
//			ImageIO.write( makeBImg(pxArrHigh), "png", new File(OUTFILE+"_high.png"));
			
//			//Here's the odd part
//			Filterer filterOdd = (W, H, i, j, clr) -> i%2==1 && j%2==1;
//			Complex[][][] pxArrOdd = FFTAll(applyFilter(wvArr, filterOdd), false);
//			makeBright(pxArrOdd, targetBrightness);
//			ImageIO.write( makeBImg(pxArrOdd), "png", new File(OUTFILE+"_Base"+System.currentTimeMillis()+".png"));
			
//			//Here's the even part
//			Filterer filterEven = (W, H, i, j, clr) -> i%2==0 && j%2==0;
//			Complex[][][] pxArrEven = FFTAll(applyFilter(wvArr, filterEven), false);
//			makeBright(pxArrEven, targetBrightness);
//			ImageIO.write( makeBImg(pxArrEven), "png", new File(OUTFILE+"_Even.png"));
			
//			System.out.println("Finished printing standard filters");
			
//			//Here's a random selection of all parts
//			for(int r=0;r<10;r++) {
//				Filterer filterRandom = (W, H, i, j, clr) -> random.nextFloat() < ((double)(i+j))/(W+H);
//				Complex[][][] pxArrRandom = FFTAll(applyFilter(wvArr, filterRandom), false);
//				makeBright(pxArrRandom, targetBrightness);
//				ImageIO.write( makeBImg(pxArrRandom), "png", new File(OUTFILE+"_randomantiweighted"+System.currentTimeMillis()+".png"));
//			}
//			System.out.println("Finished printing random(s).");
			
//			//Here's one that is colored based on the wave's number
//			Complex[][][] wvArrColor = new Complex[3][W][H];
//			for(int i=0;i<W;i++)
//				for(int j=0;j<H;j++) {
//					double zz = -Math.log(gcd(i,W)*gcd(j,H)/(double)(W*H));
//					wvArrColor[0][i][j] = wvArr[0][i][j].scale(1.5);
//					wvArrColor[1][i][j] = wvArr[0][i][j].scale(0.5);
//					wvArrColor[2][i][j] = wvArr[0][i][j].scale(Math.cbrt(zz));
//				}
//			Complex[][][] pxArrColor = FFTAll(wvArrColor, false);
//			//makeBright(pxArrColor, targetBrightness);
//			ImageIO.write( makeBImg(pxArrColor), "png", new File(OUTFILE+"_Color"+System.currentTimeMillis()+".png"));

			//Here's one that is colored based on the wave's number
			Complex[][][] wvArrColor = new Complex[3][W][H];
			for(int i=0;i<W;i++)
				for(int j=0;j<H;j++) {
					double zz = -Math.log(gcd(i,W)*gcd(j,H)/(double)(W*H));
					wvArrColor[0][i][j] = wvArr[0][i][j].scale(2.0*Math.pow(zz, 0.2));
					wvArrColor[1][i][j] = wvArr[0][i][j].scale(0.4);
					wvArrColor[2][i][j] = wvArr[0][i][j].scale(0.6*Math.pow(zz, 0.8));
				}
			Complex[][][] pxArrColor = FTAll(wvArrColor, false);
			//makeBright(pxArrColor, targetBrightness);
			ImageIO.write( makeBImg(pxArrColor), "png", new File(OUTFILE+"_Colors"+System.currentTimeMillis()+".png"));
			
//			//here's printing a graph of the frequency-space:
//			makeBright(wvArr, targetBrightness);
//			BufferedImage freqBImg = makeBImg(wvArr);
//			ImageIO.write(freqBImg, "png", new File(OUTFILE+"_graph.png"));
//			System.out.println("Finished printing 'graph'");
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		System.out.println("Ending.");
	}
	
	/**
	 * Filters a specific wave-decomposition of an image using a specified filter
	 * @param wvArr wave-decomposition of an image
	 * @param filter the filter to be applied
	 * @returna new wave-decomposition that is the filter applied to wvArr
	 */
	public static Complex[][][] applyFilter(Complex[][][] wvArr, Filterer filter) {
		Complex[][][] wvArrOut = new Complex[3][W][H];
		for(int i=0;i<W;i++)
			for(int j=0;j<H;j++)
				for(int clr=0;clr<3;clr++) {
					if(filter.filter(W, H, i, j, clr))
						wvArrOut[clr][i][j] = wvArr[clr][i][j];
					else
						wvArrOut[clr][i][j] = Complex.ZERO;
				}
		return wvArrOut;
	}
	
	/**
	 * Make a buffered image from a pixel-decomposition. Only takes the real part of each r,g,b value.
	 * Requires that H and W be the true widths and heights of the image
	 * @param pxArr image in (complex) pixel array form
	 * @return the image as a buffered image
	 */
	public static BufferedImage makeBImg(Complex[][][] pxArr) {
		BufferedImage out = new BufferedImage(W, H, BufferedImage.TYPE_INT_RGB);
		for(int i=0;i<W;i++)
			for(int j=0;j<H;j++) {
				int red = (int)(Math.round(pxArr[0][i][j].r));
				int green = (int)(Math.round(pxArr[1][i][j].r));
				int blue = (int)(Math.round(pxArr[2][i][j].r));
//				System.out.println("pixel: ("+red+", "+green+", "+blue+")");
				if(red < 0) red = 0;
				if(red > 255) red = 255;
				if(green < 0) green = 0;
				if(green > 255) green = 255;
				if(blue < 0) blue = 0;
				if(blue > 255) blue = 255;
				
				int rgb = (red<<16) + (green<<8) + (blue<<0);
				out.setRGB(i, j, rgb);
			}
		return out;
	}

	// turns a bufferedImage into a complex pixel array
	/**
	 * Create the (complex) pixel array representing a given buffered image
	 * Requires that H and W be the true widths and heights of the image
	 * Values range from 0 to 255 (real)
	 * @param bimg the buffered image to represent as a pixel array
	 * @return the pixel-array representation
	 */
	public static Complex[][][] getComplexPxlArr(BufferedImage bimg) {
		Complex[][][] output = new Complex[3][W][H];
		for(int i=0;i<W;i++)
			for(int j=0;j<H;j++) {
				int rgb = bimg.getRGB(i, j);
				int red = (rgb&0xff0000)>>>16;
				int blue = (rgb&0x00ff00)>>>8;
				int green = (rgb&0x0000ff)>>>0;
				output[0][i][j] = new Complex(red, 0);
				output[1][i][j] = new Complex(blue, 0);
				output[2][i][j] = new Complex(green, 0);
			}
		return output;
	}
	
	//TODO: look up the standard way to brighten up a picture. Hopefully, only this function will have to be changed.
	/**
	 * Finds the total brightness of a picture.
	 * Always non-negative and scales linearly with the values of pxArr.
	 * @param pxArr the (complex) pixel array to find the brightness of
	 * @return the brightness of the pixel array
	 */
	public static double getBrightness(Complex[][][] pxArr) {
		double out = 0;
		for(int i=0;i<W;i++)
			for(int j=0;j<H;j++)
				for(int clr=0;clr<3;clr++) {
					out += pxArr[clr][i][j].norm();
				}
		return out;
	}
	/**
	 * Scales the (r,g,b values of) pixels in a pixel array so that its brightness will be the target brightness
	 * Modifies the pixel array in-place.
	 * @param pxArr the (complex) pixel array to modify the brightness of
	 * @param targetBrightness the target brightness for the pixel array
	 */
	public static void makeBright(Complex[][][] pxArr, double targetBrightness) {
		double firstBrightness = getBrightness(pxArr);
		double ratio = targetBrightness / firstBrightness;
		for(int i=0;i<W;i++)
			for(int j=0;j<H;j++)
				for(int clr=0;clr<3;clr++) {
					pxArr[clr][i][j] = pxArr[clr][i][j].scale(ratio);
				}
	}
	
	
	/**
	 * Finds the gcd of two nonnegative integers using the Euclidean method
	 * @param n
	 * @param m
	 * @return the gcd of n and m
	 */
	public static int gcd(int n, int m) {
		if(n == 0)
			return m;
		else if(m == 0)
			return n;
		else if(n > m)
			return gcd(n%m, m);
		else
			return gcd(n, m%n);
	}
	
	
	//complex exponentials for all powers of a primitive Mth root ot unity
	private static Complex[] expsM;
	/**
	 * Utility function for computing a unit-norm complex number that is, roughly, a (num/dem)th root of unity
	 * @param num
	 * @param den
	 * @return exp(2*pi*i* (num/den) )
	 */
	public static Complex getExpsM(int num, int den) {
		int i = expsM.length * num / den; // this is the position within expsM
		return expsM[((i%expsM.length)+expsM.length)%expsM.length];
	}
	
	/**
	 * Performs a Discrete Fourier Transform (DFT) on a given input. This is not optimized in any way, and henceforth the Slow Fourier Transfrom (SFT)
	 * There are no restrictions on W or H while doing this.
	 * Requires an input of size W by H
	 * @param in A (complex) array on which to perform the DFT.
	 * @param forward Wwhether to do a forward DFT or backwards DFT
	 * @return A new array that is in's DFT.
	 */
	public static Complex[][] SFT(Complex[][] in, boolean forward) {
		Complex[][] out = new Complex[W][H];
		for(int i=0;i<W;i++)
			for(int j=0;j<H;j++) {
				
				// set out[i][j];
				out[i][j] = Complex.ZERO;
				for(int g=0;g<W;g++)
					for(int k=0;k<H;k++) {
						//out[i][j] += in[g][k]*exp(-ig/w)*exp(-jk/h)
						out[i][j] = Complex.add(out[i][j], 
								Complex.multiply(in[g][k],
										getExpsM((g*i*H + k*j*W) * (forward?-1:1), H*W)
										)
								);
					}
			}
		return out;
	}
	/**
	 * Performs a Discrete Fourier Transform (DFT) on a given input.
	 * This is optimized with the FFT algorithm, making it the only viable choice for large images.
	 * Requires W and H to be powers of 2
	 * @param in a (complex) array on which to perform the DFT.
	 * @param sizeX the first dimension of the array
	 * @param sizeY the second dimension of the array
	 * @param forward whether to do a forward DFT or backwards DFT
	 * @return a new array that is in's DFT.
	 */
	public static Complex[][] FFT1(Complex[][] in, int sizeX, int sizeY, boolean forward) {
		if(sizeX == 1 && sizeY == 1) {
			Complex[][] out = new Complex[1][1];
			out[0][0] = in[0][0];
			return out;
		} else if(sizeX > 1) { //halve the x-axis
			Complex[][] inEV = new Complex[sizeX/2][sizeY];
			for(int i=0;i<sizeX/2;i++)
				for(int j=0;j<sizeY;j++) {
					inEV[i][j] = in[2*i][j];
				}
			Complex[][] inOD = new Complex[sizeX/2][sizeY];
			for(int i=0;i<sizeX/2;i++)
				for(int j=0;j<sizeY;j++) {
					inOD[i][j] = in[2*i+1][j];
				}
			
			Complex[][] outEV = FFT1(inEV, sizeX/2, sizeY, forward);
			Complex[][] outOD = FFT1(inOD, sizeX/2, sizeY, forward);
			
			Complex[][] out = new Complex[sizeX][sizeY];
			for(int i=0;i<sizeX;i++)
				for(int j=0;j<sizeY;j++) {
					Complex factor = getExpsM(i * (forward?-1:1), sizeX);					
					
					//out[i][j] = outEV[i][j] + outOD[i][j] * exp(r/sizeX * forward?1:-1)
					out[i][j] = Complex.add(outEV[i%(sizeX/2)][j], Complex.multiply(outOD[i%(sizeX/2)][j], factor));
				}
			
			return out;
		} else { //halve the y-axis
			Complex[][] inEV = new Complex[sizeX][sizeY/2];
			for(int i=0;i<sizeX;i++)
				for(int j=0;j<sizeY/2;j++) {
					inEV[i][j] = in[i][2*j];
				}
			Complex[][] inOD = new Complex[sizeX][sizeY/2];
			for(int i=0;i<sizeX;i++)
				for(int j=0;j<sizeY/2;j++) {
					inOD[i][j] = in[i][2*j+1];
				}
			
			Complex[][] outEV = FFT1(inEV, sizeX, sizeY/2, forward);
			Complex[][] outOD = FFT1(inOD, sizeX, sizeY/2, forward);
			
			Complex[][] out = new Complex[sizeX][sizeY];
			for(int i=0;i<sizeX;i++)
				for(int j=0;j<sizeY;j++) {
					Complex factor = getExpsM(j * (forward?-1:1), sizeY);					
					//out[i][j] = outEV[i][j] + outOD[i][j] * exp(s/sizeY * forward?1:-1)
					out[i][j] = Complex.add(outEV[i][j%(sizeY/2)], Complex.multiply(outOD[i][j%(sizeY/2)], factor));
				}
			
			return out;
		}
	}
	/**
	 * Performs a Discrete Fourier Transform (DFT) on a given input.
	 * This is further otimized than FFT1, about twice as fast for medium images.
	 * Requires in to be a W by H array.
	 * @param in a (complex) array on which to perform the DFT.
	 * @param startX the first x index to consider
	 * @param deltaX the y step-size
	 * @param sizeX the first dimension of the array
	 * @param startY the first y index to consider
	 * @param deltaY the x step-size
	 * @param sizeY the second dimension of the array
	 * @param forward whether to do a forward DFT or backwards DFT
	 * @return a new array that is in's DFT.
	 */
	//NOTE that this handles any order of dividing coordinates, but it's easiest to divide sizeX to 1 first.
	public static Complex[][] FFT2(Complex[][] in, int startX, int deltaX, int sizeX, int startY, int deltaY, int sizeY, boolean forward) {
		Complex[][] out = new Complex[sizeX][sizeY];
		if(sizeX==1 && sizeY==1) { // base case - you don't have anything here.
			out[0][0] = in[startX][startY];
			return out;
		} else if(sizeX==1 && sizeY>1){ // in this case, divide the y-coordinate by 2:
			//get even and odd parts:
			Complex[][] outEVEN = FFT2(in, startX, deltaX, sizeX, startY, deltaY*2, sizeY/2, forward);
			Complex[][] outODD  = FFT2(in, startX, deltaX, sizeX, startY+deltaY, deltaY*2, sizeY/2, forward);
			//build the out[][]:
			for(int i=0;i<sizeX;i++)
				for(int j=0;j<sizeY;j++) {
					out[i][j] = Complex.add(outEVEN[i][j%(sizeY/2)], 
							Complex.multiply(outODD[i][j%(sizeY/2)], getExpsM(forward?j:-j, sizeY))
							);
				}
			return out;
		} else { // in this case, divide the x-coordinate by 2:
			//get even and odd parts:
			Complex[][] outEVEN = FFT2(in, startX, deltaX*2, sizeX/2, startY, deltaY, sizeY, forward);
			Complex[][] outODD  = FFT2(in, startX+deltaX, deltaX*2, sizeX/2, startY, deltaY, sizeY, forward);
			//build the out[][]:
			for(int i=0;i<sizeX;i++)
				for(int j=0;j<sizeY;j++) {
					out[i][j] = Complex.add(outEVEN[i%(sizeX/2)][j], 
							Complex.multiply(outODD[i%(sizeX/2)][j], getExpsM(forward?i:-i, sizeX))
							);
				}
			return out;
		}
	}
	/**
	 * Performs a DFT on each color individually of a pixel array.
	 * May call to SFT, FFT1, or FFT2 with little modification.
	 * A general all-purpose tool to do obtain the fourier coefficients of a pixel array
	 * @param in a (complex) array on which to perform the DFT.
	 * @param forward whether to do a forward DFT or backwards DFT
	 * @return a (complex) wave array
	 */
	public static Complex[][][] FTAll(Complex[][][] in, boolean forward) {
		Complex[][][] out = new Complex[3][][];
		for(int clr=0;clr<3;clr++) {
//			out[clr] = SFT(in[clr], forward);
//			out[clr] = FFT1(in[clr], W, H, forward);
			out[clr] = FFT2(in[clr], 0, 1, W, 0, 1, H, forward);
		}
		if(forward) { // if we are doing a forward transform, we need to scale. By definition, we don't for reverse
			for(int i=0;i<W;i++)
				for(int j=0;j<H;j++)
					for(int clr=0;clr<3;clr++) {
						out[clr][i][j] = out[clr][i][j].scale(1.0/(W*H));
					}
		}
		return out;
	}
}
