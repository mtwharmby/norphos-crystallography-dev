package uk.co.norphos.crystallography.toolkit;

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
