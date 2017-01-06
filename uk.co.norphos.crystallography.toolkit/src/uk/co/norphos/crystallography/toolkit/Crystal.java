package uk.co.norphos.crystallography.toolkit;

public class Crystal {
	
	private UnitCell unitCell;
	
	public Crystal(Lattice lattice) {
		unitCell = new UnitCell(lattice);
	}
	
	public UnitCell getUnitCell() {
		return unitCell;
	}

}
