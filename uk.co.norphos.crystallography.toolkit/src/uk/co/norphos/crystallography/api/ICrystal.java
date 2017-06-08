package uk.co.norphos.crystallography.api;

import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * FIXME
 * @author Michael Wharmby
 *
 */
public interface ICrystal {
	
	/**
	 * Return the current unit cell (including the {@link Lattice}) of this {@
	 * link ICrystal}.
	 * @return IUnitCell current unit cell of this crystal
	 */
	IUnitCell getUnitCell();
	
	/**
	 * Update the unit cell (and {@link Lattice}) of this {@link ICrystal}.
	 * @param unitCell IUnitCell describing the new lattice
	 */
	void setUnitCell(IUnitCell unitCell);
	
	/**
	 * Calculate and return a list of the diffracted beams 
	 * ({@link MillerPlane}s) which will be observed up to a d-spacing limit.
	 * @param dSpacing Double limit up to which reflections could be observed. 
	 * If null, returns current list.
	 * @return List<MillerPlane> expected reflections up to d-space limit
	 * TODO should be a Map?
	 */
	List<MillerPlane> calculateExpectedReflections(Double dSpacing);
	
	/**
	 * Return the current list of diffracted beams which will be observed up 
	 * to the previously set d-spacing limit. If no limit has been previously 
	 * set (i.e. {@link #calculateExpectedReflections(Double)} has not been 
	 * called), an empty list may be returned.
	 * @return List<MillerPlane> expected reflections
	 * TODO should be a Map?
	 */
	default List<MillerPlane> getExpectedReflections() {
		return calculateExpectedReflections(null);
	}
	
	/**
	 * Return the current {@link ISpaceGroup} of this crystal.
	 * @return {@link ISpaceGroup}
	 */
	ISpaceGroup getSpaceGroup();
	
	/**
	 * Change the symmetry of this crystal by updating the {@link ISpaceGroup}.
	 * @param spaceGroup {@link ISpaceGroup}
	 */
	void setSpaceGroup(ISpaceGroup spaceGroup);
	
	/**
	 * Return the list of symmetry unique atoms which make up this crystal 
	 * structure.
	 * @return List<Atom> list of atoms forming the asymmetric unit
	 */
	List<Atom> getAsymmetricUnit();
	
	void setAsymmetricUnit(List<Atom> atoms);
	
	/**
	 * Add atom to the asymmetric unit.
	 * @param atom {@link Atom}
	 */
	default void addAtom(Atom atom) {
		getAsymmetricUnit().add(atom);
	}
	
	/**
	 * Remove an atom from the asymmetric unit.
	 * @param atom {@link Atom}
	 */
	default void removeAtom(Atom atom) {
		getAsymmetricUnit().remove(atom);
	}
	
	/**
	 * Remove an atom from the asymmetric unit by its index.
	 * @param i int index of {@link Atom} to remove
	 */
	default void removeAtom(int i) {
		getAsymmetricUnit().remove(i);
	}
	
	/**
	 * Return the UB matrix of this {@link ICrystal}.
	 * @return RealMatrix
	 */
	RealMatrix getOrientationMatrix();
	
	/**
	 * Change the UB matrix of this {@link ICrystal}.
	 * @param ubMatrix RealMatrix
	 */
	void setOrientationMatrix(RealMatrix ubMatrix);
	
	/**
	 * Calculate the density of the crystal structure, determined from the 
	 * mass of the contents of the unit cell and its volume.
	 * @return double density in g cm<sup>-3</sup>
	 */
	double calculateDensity();

}
