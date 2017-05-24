package uk.co.norphos.crystallography.api;

/**
 * A generic coordinate in up to three dimensions.
 * 
 * @author Michael Wharmby
 *
 */
public interface ICoord {
	
	/**
	 * Return the first dimension's magnitude 
	 */
	double getX();
	
	/**
	 * Return the second dimension's magnitude 
	 */
	double getY();
	
	/**
	 * Return the third dimension's magnitude 
	 */
	double getZ();
	

}
