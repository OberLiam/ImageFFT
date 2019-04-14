package main;

/**
 * @author Liam
 * Copmlex Numbers
 */
public class Complex {
	public static final Complex ZERO = new Complex(0,0);
	public static final Complex ONE = new Complex(1,0);
	
	// the real and imaginary parts to a function
	public double r,i;
	
	public Complex(double r, double i) {
		this.r = r;
		this.i = i;
	}
	
	/**
	 * Returns a complex number with the given magnitude and phase
	 * @param R magnitude
	 * @param theta phase/complex argument
	 * @return the specified complex number
	 */
	public static Complex makeRT(double R, double theta) {
		return new Complex(R*Math.cos(theta), R*Math.sin(theta));
	}
	
	/**
	 * Multiply two complex numbers, a and b
	 * @param a
	 * @param b
	 * @return a*b
	 */
	public static Complex multiply(Complex a, Complex b) {
		return new Complex(a.r*b.r - a.i*b.i, a.r*b.i + a.i*b.r);
	}
	/**
	 * Multiply two complex numbers, a and b
	 * @param a
	 * @param b
	 * @return a*b
	 */
	public static Complex add(Complex a, Complex b) {
		return new Complex(a.r + b.r, a.i + b.i);
	}
	public static Complex addAll(Complex... as) {
		Complex out = Complex.ZERO;
		for(Complex a : as)
			out = add(out, a);
		return out;
	}
	
	public Complex scale(double x) {
		return new Complex(x*r, x*i);
	}
	public Complex con() {
		return new Complex(r,-i);
	}
	public double norm() {
		return Math.hypot(r, i);
	}
	public double arg() {
		return Math.atan2(i, r);
	}
}
