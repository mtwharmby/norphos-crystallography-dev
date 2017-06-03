package uk.co.norphos.crystallography.api;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * IUnitCell describes the size of the repeating 3d tile of a crystal. It 
 * consists of a {@link Lattice} and provides methods to change the lattice 
 * parameters. Furthermore it allows the calculation of values derived or 
 * dependent on from the lattice (e.g. volume or d-spacing). 
 * 
 * @author Michael Wharmby
 *
 */
public interface IUnitCell extends Comparable<IUnitCell> {
	
	/**
	 * Return the real-space lattice parameters for this IUnitCell.
	 * @return {@link Lattice}
	 */
	Lattice getLattice();
	
	/**
	 * Change the real-space {@link Lattice} of this IUnitCell. 
	 * 
	 * @param {@link Lattice} lattice
	 */
	void setLattice(Lattice lattice);
	
	/**
	 * Set one or more of the lattice parameters. A new lattice instance will 
	 * be created and all values will be recalculated.
	 * 
	 * @param a Double length in Angstroms
	 * @param b Double length in Angstroms
	 * @param c Double length in Angstroms
	 * @param alpha Double angle in degrees
	 * @param beta Double angle in degrees
	 * @param gamma Double angle in degrees
	 * @throws LatticeException if the the given lattice parameters are 
	 *         incompatible or inconsistent with the {@link CrystalSystem} 
	 */
	void setLattice(Double a, Double b, Double c, Double alpha, Double beta, Double gamma) throws LatticeException;
	
	/**
	 * Return the reciprocal-space lattice parameters for this IUnitCell.
	 * @return {@link Lattice}
	 */
	default Lattice getReciprocalLattice() {
		return getReciprocal().getLattice();
	}
	
	/**
	 * Return the {@link CrystalSystem} (i.e. lattice system) of the real-
	 * space lattice. N.B. This may differ from the {@link CrystalSystem} 
	 * reported from symmetry, e.g. the hexagonal lattice system contains 
	 * both hexagonal and trigonal lattices; the rhombohedral lattice 
	 * system is incorporated into the trigonal crystal system.
	 * @return {@link CrystalSystem}
	 */
	CrystalSystem getLatticeSystem();
	
	/**
	 * Return the volume of the unit cell.
	 * @return double volume in Angstrom^3
	 */
	double getVolume();
	
	/**
	 * Return the metric tensor (G-matrix) for the real-space unit cell.
	 * @return RealMatrix G-matrix
	 */
	RealMatrix getMetricTensor();
	
	/**
	 * Return the metric tensor of the reciprocal-space unit cell.
	 * @return RealMatrix reciprocal-space G-matrix
	 */
	default RealMatrix getReciprocalMetricTensor() {
		return getReciprocal().getMetricTensor();
	}
	
	/**
	 * Return the reciprocal-space equivalent of this IUnitCell.
	 * @return IUnitCell
	 */
	IUnitCell getReciprocal();
	
//FIXME
//	/**
//	 * Convert a vector in Cartesian coordinates to its equivalent in the 
//	 * fractional coordinate system of this unit cell.
//	 * 
//	 * @param cartVector Vector3D in Cartesian coordinates
//	 * @return Vector3D in fractional coordinates of the current lattice
//	 */
//	default Vector3D fractionalize(Vector3D cartVector) {
//		//TODO
//		return null;
//	}
//	
//	/**
//	 * Convert a vector in fractional coordinates of this unit cell into an 
//	 * equivalent vector in Cartesian coordinates.
//	 * 
//	 * @param fracVector Vector3D in fractional coordinates
//	 * @return Vector3D in Cartesian coordinates
//	 */
//	default Vector3D orthogonalize(Vector3D fracVector) {
//		//TODO
//		return null;
//	}
//	
//	/**
//	 * Return matrix to convert Cartesian coordinates into fractional 
//	 * coordinates for this unit cell's lattice.
//	 * @return RealMatrix
//	 */
//	RealMatrix getFractionalizationMatrix();
//	
//	/**
//	 * Return matrix to convert fractional coordinates of this unit cell's 
//	 * lattice into Cartesian coordinates.
//	 * @return RealMatrix
//	 */
//	RealMatrix getOrthogonalizationMatrix();
//	
//	/**
//	 * Calculate the length of a vector specified in fractional coordinates of 
//	 * this unit cell.
//	 * 
//	 * @param fracVec Vector3D in fractional coordinates
//	 * @return double length of vector
//	 */
//	double calculateLength(Vector3D fracVec);
//	
//	/**
//	 * Calculate the distance between two sites specified in fractional 
//	 * coordinates of this unit cell.
//	 * 
//	 * @param site1 Vector3D in fractional coordinates
//	 * @param site2 Vector3D in fractional coordinates
//	 * @return
//	 */
//	default double calculateDistance(Vector3D site1, Vector3D site2) {
//		//TODO
//		return 0;
//	}
//	
//	/**
//	 * Calculate the angle between two vectors specified in fractional 
//	 * coordinates of this unit cell.
//	 * 
//	 * @param fracVec1 Vector3D in fractional coordinates
//	 * @param fracVec2 Vector3D in fractional coordinates
//	 * @return double angle between vectors in radians
//	 */
//	double calculateAngle(Vector3D fracVec1, Vector3D fracVec2);
//	
//	/**
//	 * Calculate the angle between site 1 and site3 at site2 (i.e. the angle 
//	 * between the vectors site1-site2 and site2-site3, c.f. bond angle).
//	 * 
//	 * @param site1 Vector3D in fractional coordinates
//	 * @param site2 Vector3D in fractional coordinates
//	 * @param site3 Vector3D in fractional coordinates
//	 * @return double angle at site2 in radians
//	 */
//	default double calculateAngle(Vector3D site1, Vector3D site2, Vector3D site3) {
//		//TODO
//		return 0;
//	}
//	
//	/**
//	 * Calculate the angle between the planes containing site1, site2 and 
//	 * site3 and site2, site3 and site4.
//	 * 
//	 * @param site1 Vector3D in fractional coordinates
//	 * @param site2 Vector3D in fractional coordinates
//	 * @param site3 Vector3D in fractional coordinates
//	 * @param site4 Vector3D in fractional coordinates
//	 * @return double angle between planes in radians
//	 */
//	default double calculateDihedralAngle(Vector3D site1, Vector3D site2, Vector3D site3, Vector3D site4) {
//		//TODO
//		return 0;
//	}
//	
//	/**
//	 * Maximum {@link MillerIndex} for the given d-spacing limit.
//	 * 
//	 * @param dSpacing in Angstrom^-1 TODO Right?
//	 * @return MillerIndex maximum hkl observable
//	 */
//	MillerIndex getMaxMillerIndex(double dSpacing);
//	
//	/**
//	 * Return the d-space value for a specific {@link MillerIndex}.
//	 *  
//	 * @param hkl {@link MillerIndex}
//	 * @return double d-spacing in Angstrom^-1 TODO Right?
//	 */
//	double getDSpacing(MillerIndex hkl);
//	
//	/**
//	 * Determine whether this IUnitCell is similar to another one, within certain tolerances. 
//	 *  
//	 * @param other IUnitCell to compare
//	 * @param lengthTol Double length tolerance (if null, default to 0.02 - 2%)
//	 * @param angleTol Double angle tolerance (if null, default to 1degree)
//	 * @return boolean true if this and other are same within tolerance
//	 */
//	boolean isSimilar(IUnitCell other, Double lengthTol, Double angleTol);
}
