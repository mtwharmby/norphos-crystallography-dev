package uk.co.norphos.crystallography.api;

/**
 * IUnitCell describes the size of the repeating 3d tile of a crystal. It 
 * consists of a {@link Lattice} and provides methods to change the lattice 
 * parameters. Furthermore it allows the calculation of values derived or 
 * dependent on from the lattice (e.g. volume or d-spacing). 
 * 
 * @author Michael Wharmby
 * 
 * TODO Should/could we specify units for our lattice params? JScience?
 *
 */
public interface IUnitCell {
	
	Lattice getLattice();
	
	void setLattice(Lattice lattice);
	
	/**
	 * Set one or more of the lattice parameters. If the current lattice should 
	 * be retained and updated, set update true. Otherwise the 
	 * {@link CrystalSystem} of the IUnitCell should be redetermined.
	 * @param a Double a length in angstroms
	 * @param b Double b length in angstroms
	 * @param c Double c length in angstroms
	 * @param alpha Double alpha angle in degrees
	 * @param beta Double beta angle in degrees
	 * @param gamma Double gamma angle in degrees
	 * @param update boolean, if true lattice parameters will be updated; if 
	 *        false a new {@link CrystalSystem} will be redetermined
	 * @throws LatticeException if the the given lattice parameters are 
	 *         incompatible with the current {@link CrystalSystem} 
	 */
	void setLatticeParameters(Double a, Double b, Double c, Double alpha, Double beta, Double gamma, boolean update) throws LatticeException;
	
	/**
	 * Return the a lattice parameter in angstroms
	 */
	double getA();
	
	/**
	 * Set the a lattice parameter (in angstroms)
	 */
	void setA(double a);
	
	/**
	 * Return the value of lattice a parameter
	 */
	double getB();
	
	/**
	 * Set value of lattice b parameter
	 */
	void setB(double b);
	
	/**
	 * Return the value of lattice c parameter
	 */
	double getC();
	
	/**
	 * Set value of lattice c parameter
	 */
	void setC(double c);
	
	/**
	 * Return the alpha lattice parameter (in degrees) 
	 */
	double getAlpha();
	
	/**
	 * Set value of lattice alpha parameter in degrees
	 */
	void setAlpha(double alpha);
	
	/**
	 * Return the beta lattice parameter (in degrees)
	 */
	double getBeta();
	
	/**
	 * Set value of lattice beta parameter in degrees
	 */
	void setBeta(double beta);
	
	/**
	 * Return the gamma lattice parameter (in degrees)
	 */
	double getGamma();
	
	/**
	 * Set value of lattice gamma parameter in degrees
	 */
	void setGamma(double gamma);
	
	/**
	 * Return the volume of the current lattice 
	 */
	double getVolume();
	
	/**
	 * Return the {@link CrystalSystem} determined for the current lattice
	 */
	CrystalSystem getCrystalSystem();
	
	/**
	 * Calculate the distance between any two points in the lattice of this 
	 * {@link IUnitCell}, specified in fractional coordinates of the lattice.
	 * 
	 * @param point1 {@link ICoord} first point in fractional coordinates
	 * @param point2 {@link ICoord} second point in fractional coordinates
	 * @return double distance in angstroms
	 */
	double calculateDistanceBetweenTwoPoints(ICoord point1, ICoord point2);
	
	/**
	 * Calculate the distance in angstroms between two Miller planes with the 
	 * same hkl indices
	 * 
	 * @param hkl {@link ICoord} containing hkl indices
	 * @return double distance in angstrom between Miller planes
	 */
	double calculateDSpacing(ICoord hkl);
	
}
