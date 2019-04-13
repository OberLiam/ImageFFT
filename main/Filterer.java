package main;

/**
 * @author Liam
 * This interface allows for arbitrary (wave-based) filters to be a applied to an image
 */
public interface Filterer {
	public boolean filter(int w, int h, int i, int j, int clr);
}
