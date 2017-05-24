package uk.co.norphos.crystallography.toolkit;

import uk.co.norphos.crystallography.api.CrystalSystem;
import uk.co.norphos.crystallography.api.Lattice;

public class Crystal {
	
	private UnitCell unitCell;
	private CrystalSystem crystalSystem;
	
	public Crystal(Lattice lattice) {
		unitCell = new UnitCell(lattice);
	}
	
	public UnitCell getUnitCell() {
		return unitCell;
	}

}
