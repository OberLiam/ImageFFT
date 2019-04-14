package main;

/**
 * @author Liam
 * This interface specifies an arbitrary (wave-based) filter that can be a applied to an image
 */
public interface Filterer {
	/**
	 * @param w width of the image
	 * @param h height of the image
	 * @param i x-coordinate of the wave
	 * @param j y-coordinate of the wave
	 * @param clr color of the wave (0 = red, 1 = green, 2 = blue)
	 * @return whether the the specified wave should be kept (i.e. true => keep, false => discard)
	 */
	public boolean filter(int w, int h, int i, int j, int clr);
}
